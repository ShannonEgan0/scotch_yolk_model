import tkinter as tk

import matplotlib.pyplot
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation


# Main function
def main():
    root = tk.Tk()
    root.title('Piston Simulation')
    MainWindow(root).grid(row=0, column=0)
    root.mainloop()


# Differentiation function, for velocity and acceleration calculation
def inst_v(data, time):
    y1 = data[-2]
    y2 = data[-1]
    t1 = time[-2]
    t2 = time[-1]
    v = (y2 - y1) / (t2 - t1)
    return v


# Function to determine crank location for straight shaft
def piston(angular_velocity, radius, time, crank_length=1, initial_angle=0, yolk_angle=0, scotch=False):
    angle = initial_angle + time * angular_velocity     # Angle in degrees from initial_angle at time
    flywheel_position = [radius * np.sin(angle), radius * np.cos(angle)]  # Connection cartesian position
    if scotch:
        piston_position = flywheel_position[0] + flywheel_position[1] * np.tan(np.deg2rad(yolk_angle))
    else:
        piston_position = flywheel_position[0] + np.sqrt(crank_length ** 2 - flywheel_position[1] ** 2) - crank_length
    return angle, flywheel_position, piston_position


def ang_velocity(freq):
    w = 2 * np.pi / (1 / freq)
    stept = 0.01 * (2 * np.pi / (1 / 0.5)) / (2 * np.pi / (1 / freq))
    return w, stept


class MainWindow(tk.Frame):
    def __init__(self, parent, running=False):
        tk.Frame.__init__(self, parent)

        self.canvas = None
        self.canvast = None
        self.anit = None
        self.ani = None
        self.running = running
        self.parent = parent

        self.points = 3000

        self.paramframe = tk.Frame()
        self.paramframe.grid(row=0, column=0, sticky=tk.W)

        self.plotframe = tk.Frame()
        self.plotframe.grid(row=1, column=0)

        self.rad_label = tk.Label(self.paramframe, text="Crank Radius (m): ")
        self.rad_label.grid(row=0, column=0, sticky=tk.E)
        self.rad_entry = tk.Entry(self.paramframe)
        self.rad_entry.insert(0, "1")
        self.rad_entry.grid(row=0, column=1)

        self.freq_label = tk.Label(self.paramframe, text="Motor Frequency (Hz): ")
        self.freq_label.grid(row=1, column=0, sticky=tk.E)
        self.freq_entry = tk.Entry(self.paramframe)
        self.freq_entry.insert(0, "0.5")
        self.freq_entry.grid(row=1, column=1)

        self.sc_ang_label = tk.Label(self.paramframe, text="Yoke Angle (Degrees): ")
        self.sc_ang_label.grid(row=3, column=0, sticky=tk.E)
        self.sc_ang_entry = tk.Entry(self.paramframe)
        self.sc_ang_entry.insert(0, "0")
        self.sc_ang_entry.grid(row=3, column=1)

        self.cr_len_label = tk.Label(self.paramframe, text="Connecting Rod Length (m): ")
        self.cr_len_label.grid(row=2, column=0, sticky=tk.E)
        self.cr_len_entry = tk.Entry(self.paramframe)
        self.cr_len_entry.insert(0, "10")
        self.cr_len_entry.grid(row=2, column=1)

        self.sc_button = tk.Button(self.paramframe, text="Start Animation",
                                   command=lambda: self.start_animation(float(self.rad_entry.get()),
                                                                        float(self.freq_entry.get()),
                                                                        float(self.sc_ang_entry.get()),
                                                                        float(self.cr_len_entry.get())))
        self.sc_button.grid(row=4, column=0)

        self.stop_button = tk.Button(self.paramframe, text="Stop", command=lambda: self.stop_())
        self.stop_button.grid(row=4, column=1)

        self.Rvar = tk.StringVar(value="SY")
        self.RSCrank = tk.Radiobutton(self.paramframe, text="Straight Connection", variable=self.Rvar, value="SC")
        self.RScotch = tk.Radiobutton(self.paramframe, text="Scotch Yoke", variable=self.Rvar, value="SY")
        self.RSCrank.grid(row=0, column=2, sticky=tk.W)
        self.RScotch.grid(row=1, column=2, sticky=tk.W)

        self.fig = Figure()
        self.figt = Figure()

    def stop_(self):
        if self.running:
            self.ani.pause()
            self.anit.pause()
            self.figt.clf()
            self.fig.clf()
            self.running = False

    def start_animation(self, r, f, ya, cr):
        self.stop_()
        state = self.Rvar.get()
        if state == "SY":
            self.start_scotch(r, f, ya)
        elif state == "SC":
            self.start_straight(r, f, cr)

    def start_straight(self, radius, freq, crank_length, comparison_sine=True):
        self.running = True

        # Angular velocity calculation based on frequency input
        w, stept = ang_velocity(freq)

        # Creation of a lot of empty lists, surely these can be minimized
        lengths = []
        ts, ps, ps2 = [], [], []
        first_loop = [[], []]
        vts, vps, aps, ats = [], [], [], []

        # For comparison sine wave
        cvps = []
        caps = []

        for x in range(self.points):
            t = x * stept
            a, pp, length = piston(w, radius, t, crank_length=crank_length, scotch=False)
            ts.append(t)
            ps.append(pp)
            ps2 = list(zip(*ps))
            lengths.append(length)

            # Separately recording first loop so that the limits on the scale can be inferred
            # Eliminates the requirement for testing for new limits throughout entire model
            if a <= 360:
                first_loop[0].append(pp[0])
                first_loop[1].append(pp[1])

            # Waits until sufficient recordings made for differentiation for velocity
            if len(ts) >= 2:
                vps.append(inst_v(lengths, ts))
                vts.append(t)

            # Waits until sufficient recordings are made for differentiation for acceleration
            if len(vps) >= 2:
                aps.append(inst_v(vps, ts))
                ats.append(t)

            if comparison_sine:
                if len(ts) >= 2:
                    cvps.append(inst_v(ps2[0], ts))
                if len(cvps) >= 2:
                    caps.append(inst_v(cvps, ts))

        ax = self.fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plotframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)

        at = self.figt.add_subplot()
        self.canvast = FigureCanvasTkAgg(self.figt, master=self.plotframe)
        self.canvast.get_tk_widget().grid(row=0, column=2)

        ax.axis((-radius - 0.1, radius + crank_length + 0.1, -radius - 0.1, radius + 0.1))

        mot_rot, = ax.plot([], [], color='b')
        shaft_ln, = ax.plot([], [], color='r')
        mot_rod = ax.scatter([], [], color='b')
        velocity, = at.plot([], [], color='r', label="Velocity (m/s)")
        acceleration, = at.plot([], [], color='g', label="Acceleration (m/s²)")
        displacement, = at.plot([], [], color='b', label="Displacement (m)")
        sine_d, = at.plot([], [], color='grey', label="Perfect Sine")
        sine_v, = at.plot([], [], color='grey')
        sine_a, = at.plot([], [], color='grey')
        handles, labels = at.get_legend_handles_labels()
        at.legend(handles, labels, loc='lower right')
        at.grid(which='Both')

        # Attempt to improve performance only
        ps = np.array(ps)
        lengths = np.array(lengths)
        vts = np.array(vts)
        ats = np.array(ats)

        if comparison_sine:
            atymax = max(max(aps), max(vps), max(lengths), max(caps), max(cvps), max(ps2[0])) + 0.1
        else:
            atymax = max(max(aps), max(vps), max(lengths), max(ps2[0])) + 0.1
        atymin = min(min(aps), min(vps), min(lengths)) - 0.1
        at.set_ylim(atymin, atymax)
        self.figt.tight_layout()
        self.fig.tight_layout()

        def update(frame):
            frame = int(frame)
            at.set_xlim(0, ts[frame])
            if frame < len(first_loop[0]):
                mot_rot.set_data(first_loop[0][:frame], first_loop[1][:frame])
            shaft_ln.set_data([ps[frame][0], lengths[frame] + crank_length], [ps[frame][1], 0])
            mot_rod.set_offsets((ps[frame][0], ps[frame][1]))

            velocity.set_data(vts[:frame], vps[:frame])
            acceleration.set_data(ats[:frame], aps[:frame])
            displacement.set_data(ts[:frame], lengths[:frame])

            if comparison_sine:
                sine_d.set_data(ts[:frame], ps2[0][:frame])
                sine_a.set_data(ats[:frame], caps[:frame])
                sine_v.set_data(vts[:frame], cvps[:frame])

        self.ani = FuncAnimation(self.fig, update, frames=self.points-4, interval=10)
        self.anit = FuncAnimation(self.figt, update, frames=self.points-4, interval=10)
        self.canvas.draw()
        self.canvast.draw()

    def start_scotch(self, radius, freq, ya):
        self.running = True

        fig = Figure()
        ax = fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plotframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)

        figt = Figure()
        at = figt.add_subplot()
        at.grid(which='Both')
        self.canvast = FigureCanvasTkAgg(figt, master=self.plotframe)
        self.canvast.get_tk_widget().grid(row=0, column=2)

        fig.tight_layout()
        figt.tight_layout()

        # Angular velocity and time interval calculation
        w, stept = ang_velocity(freq)

        length_with_yolk = []
        ts, ps = [], []
        first_loop = [[], []]

        vts, vps, aps, ats = [], [], [], []

        for x in range(self.points):
            t = x * stept
            a, pp, length = piston(w, radius, t, yolk_angle=ya, scotch=True)

            ts.append(t)
            ps.append(pp)

            extra = pp[1] * np.tan(ya)
            length = pp[0] + extra
            length_with_yolk.append(length)

            if a <= 360:
                first_loop[0].append(pp[0])
                first_loop[1].append(pp[1])

            if len(ts) >= 3:
                vps.append(inst_v(length_with_yolk, ts))
                vts.append(ts[-2])

            if len(vps) >= 3:
                aps.append(inst_v(vps, ts))
                ats.append(vts[-2])

        mot_rot, = ax.plot([], [], color='b')
        shaft_ln, = ax.plot([], [], color='r')
        mot_rod = ax.scatter([], [], color='b')
        velocity, = at.plot([], [], color='r', label="Velocity (m/s)")
        acceleration, = at.plot([], [], color='g', label="Acceleration (m/s²)")
        displacement, = at.plot([], [], color='b', label="Displacement (m)")
        handles, labels = at.get_legend_handles_labels()
        at.legend(handles, labels, loc='lower right')

        atymax = max(max(aps), max(vps), max(length_with_yolk)) + 0.1
        atymin = min(min(aps), min(vps), min(length_with_yolk)) - 0.1
        at.set_ylim(atymin, atymax)
        ax.axis((-radius - radius * np.tan(ya) - 0.1, radius + radius * np.tan(ya) + 0.1,
                 -1.5 * radius, 1.5 * radius))
        figt.tight_layout()

        def update(frame):
            frame = int(frame)
            at.set_xlim(0, ts[frame])
            if frame < len(first_loop[0]):
                mot_rot.set_data(first_loop[0][:frame], first_loop[1][:frame])
            shaft_ln.set_data([length_with_yolk[frame] + radius * np.tan(ya), length_with_yolk[frame]
                               - radius * np.tan(ya)], [-radius, radius])
            mot_rod.set_offsets((ps[frame][0], ps[frame][1]))
            if frame >= 3:
                velocity.set_data(vts[:frame], vps[:frame])
            if frame >= 4:
                acceleration.set_data(ats[:frame], aps[:frame])
            displacement.set_data(ts[:frame], length_with_yolk[:frame])

        self.ani = FuncAnimation(fig, update, frames=self.points-4, interval=10)
        self.anit = FuncAnimation(figt, update, frames=self.points-4, interval=10)
        self.canvas.draw()
        self.canvast.draw()


if __name__ == "__main__":
    main()
