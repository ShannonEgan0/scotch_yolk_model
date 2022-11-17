import tkinter as tk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation


def inst_v(data, time):
    y1 = data[-3]
    y2 = data[-1]
    t1 = time[-3]
    t2 = time[-1]
    v = (y2 - y1) / (t2 - t1)
    return v


class MainWindow(tk.Frame):
    def __init__(self, parent, Running=True):
        tk.Frame.__init__(self, parent)

        self.Running = Running
        self.parent = parent

        self.paramframe = tk.Frame()
        self.paramframe.grid(row=0, column=0, sticky=tk.W)

        self.plotframe = tk.Frame()
        self.plotframe.grid(row=1, column=0)

        self.rad_label = tk.Label(self.paramframe, text="Crank Radius (mm): ")
        self.rad_label.grid(row=0, column=0, sticky=tk.E)
        self.rad_entry = tk.Entry(self.paramframe)
        self.rad_entry.insert(0, "1")
        self.rad_entry.grid(row=0, column=1)

        self.freq_label = tk.Label(self.paramframe, text="Motor Frequency (Hz): ")
        self.freq_label.grid(row=1, column=0, sticky=tk.E)
        self.freq_entry = tk.Entry(self.paramframe)
        self.freq_entry.insert(0, "0.5")
        self.freq_entry.grid(row=1, column=1)

        self.cr_len_label = tk.Label(self.paramframe, text="Connecting Rod Length (mm): ")
        self.cr_len_label.grid(row=2, column=0, sticky=tk.E)
        self.cr_len_entry = tk.Entry(self.paramframe)
        self.cr_len_entry.insert(0, "10")
        self.cr_len_entry.grid(row=2, column=1)

        self.sc_ang_label = tk.Label(self.paramframe, text="Yoke Angle (Degrees): ")
        self.sc_ang_label.grid(row=3, column=0, sticky=tk.E)
        self.sc_ang_entry = tk.Entry(self.paramframe)
        self.sc_ang_entry.insert(0, "0")
        self.sc_ang_entry.grid(row=3, column=1)

        self.sc_button = tk.Button(self.paramframe, text="Start Animation",
                                   command=lambda: self.Start_Animation(float(self.rad_entry.get()),
                                                                        float(self.freq_entry.get()),
                                                                        float(self.sc_ang_entry.get()),
                                                                        float(self.cr_len_entry.get())))
        self.sc_button.grid(row=4, column=0)

        self.stop_button = tk.Button(self.paramframe, text="Stop", command=lambda: self.Stop_())
        self.stop_button.grid(row=4, column=1)

        self.Rvar = tk.StringVar(value="SY")
        self.RSCrank = tk.Radiobutton(self.paramframe, text="Straight Connection", variable=self.Rvar, value="SC")
        self.RScotch = tk.Radiobutton(self.paramframe, text="Scotch Yoke", variable=self.Rvar, value="SY")
        self.RSCrank.grid(row=0, column=2, sticky=tk.W)
        self.RScotch.grid(row=1, column=2, sticky=tk.W)

    def Stop_(self):
        self.Running = False
        self.ani.pause()
        self.anit.pause()

    def Start_Animation(self, r, f, ya, cr, t=0):
        state = self.Rvar.get()
        if state == "SY":
            self.Start_Scotch(r, f, ya)
        if state == "SC":
            self.Start_Straight(r, f, ya, cr)

    def Start_Straight(self, r, f, ya, cr, t=0, stept=0.01, p_sine=True):
        self.Running = True

        fig = Figure()
        ax = fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plotframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)

        figt = Figure()
        at = figt.add_subplot()
        self.canvast = FigureCanvasTkAgg(figt, master=self.plotframe)
        self.canvast.get_tk_widget().grid(row=0, column=2)

        ya = np.deg2rad(ya)

        a0 = 0
        w = 2 * np.pi / (1 / f)
        a = a0 + t * w
        a = np.deg2rad(a)
        print(str(r))
        # ly = 2*r/np.cos(a)
        lengths = []
        pp = [r * np.sin(a), r * np.cos(a)]
        ts, ps = [], []
        first_loop = [[], []]

        vts, vps, aps, ats = [], [], [], []

        # For comparison sine wave
        cvps = []
        caps = []

        for x in range(2500):

            t = t + stept
            a = a0 + t * w
            pp = [r * np.sin(a), r * np.cos(a)]
            ts.append(t)
            ps.append(pp)
            ps2 = list(zip(*ps))

            extra = np.sqrt(cr ** 2 - pp[1] ** 2)
            length = pp[0] + extra
            lengths.append(length - cr)

            if a <= 360:
                first_loop[0].append(pp[0])
                first_loop[1].append(pp[1])

            if len(ts) >= 3:
                vps.append(inst_v(lengths, ts))
                vts.append(ts[-2])

            if len(vps) >= 3:
                aps.append(inst_v(vps, ts))
                ats.append(vts[-2])

            if p_sine == True:

                if len(ts) >= 3:
                    cvps.append(inst_v(ps2[0], ts))

                if len(cvps) >= 3:
                    caps.append(inst_v(cvps, ts))

        ax.axis((-r - 0.1, r + cr + 0.1, -r - 0.1, r + 0.1))
        figt.tight_layout()
        fig.tight_layout()

        mot_rot, = ax.plot([], [], color='b')
        shaft_ln, = ax.plot([], [], color='r')
        mot_rod = ax.scatter([], [], color='b')
        velocity, = at.plot([], [], color='r')
        acceleration, = at.plot([], [], color='g')
        displacement, = at.plot([], [], color='b')
        sine_d, = at.plot([], [], color='grey')
        sine_v, = at.plot([], [], color='grey')
        sine_a, = at.plot([], [], color='grey')

        atymax = max(max(aps), max(vps), max(lengths), max(caps), max(cvps), max(ps2[0])) + 0.1
        atymin = min(min(aps), min(vps), min(lengths)) - 0.1
        at.set_ylim(atymin, atymax)

        def update(frame):
            frame = int(frame)
            at.set_xlim(0, ts[frame])
            if frame < len(first_loop[0]):
                mot_rot.set_data(first_loop[0][:frame], first_loop[1][:frame])
            shaft_ln.set_data([ps[frame][0], lengths[frame] + cr], [ps[frame][1], 0])
            mot_rod.set_offsets((ps[frame][0], ps[frame][1]))
            if frame >= 3:
                velocity.set_data(vts[:frame], vps[:frame])
                sine_v.set_data(vts[:frame], cvps[:frame])
            if frame >= 4:
                acceleration.set_data(ats[:frame], aps[:frame])
                sine_a.set_data(ats[:frame], caps[:frame])
            displacement.set_data(ts[:frame], lengths[:frame])
            sine_d.set_data(ts[:frame], ps2[0][:frame])

        self.ani = FuncAnimation(fig, update, frames=np.linspace(1, 2499, 2499), interval=50)
        self.anit = FuncAnimation(figt, update, frames=np.linspace(1, 2499, 2499), interval=50)
        self.canvas.draw()
        self.canvast.draw()


    def Start_Scotch(self, r, f, ya, t=0, stept=0.005):
        self.Running = True

        fig = Figure()
        ax = fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plotframe)
        self.canvas.get_tk_widget().grid(row=0, column=0)

        figt = Figure()
        at = figt.add_subplot()
        at.grid(which='Both')
        at.set_yticks([0, 5, 10, 15, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31])
        self.canvast = FigureCanvasTkAgg(figt, master=self.plotframe)
        self.canvast.get_tk_widget().grid(row=0, column=2)

        fig.tight_layout()
        figt.tight_layout()

        ya = np.deg2rad(ya)
        yl = r * np.tan(ya)

        a0 = 0
        w = 2 * np.pi / (1 / f)
        a = a0 + t * w
        a = np.deg2rad(a)
        print(str(r))
        # ly = 2*r/np.cos(a)
        length_with_yolk = []
        pp = [r * np.sin(a), r * np.cos(a)]
        ts, ps = [], []
        first_loop = [[], []]

        vts, vps, aps, ats = [], [], [], []

        for x in range(2500):

            t = t + stept
            a = a0 + t * w
            pp = [r * np.sin(a), r * np.cos(a)]
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
        velocity, = at.plot([], [], color='r')
        acceleration, = at.plot([], [], color='g')
        displacement, = at.plot([], [], color='b')
        sine_d, = at.plot([], [], color='grey')
        sine_v, = at.plot([], [], color='grey')
        sine_a, = at.plot([], [], color='grey')

        atymax = max(max(aps), max(vps), max(length_with_yolk)) + 0.1
        atymin = min(min(aps), min(vps), min(length_with_yolk)) - 0.1
        at.set_ylim(atymin, atymax)
        ax.axis((-r - yl - 0.1, r + yl + 0.1, -1.5 * r, 1.5 * r))

        def update(frame):
            frame = int(frame)
            at.set_xlim(0, ts[frame])
            if frame < len(first_loop[0]):
                mot_rot.set_data(first_loop[0][:frame], first_loop[1][:frame])
            shaft_ln.set_data([length_with_yolk[frame] + r * np.tan(ya), length_with_yolk[frame] - r * np.tan(ya)],
                              [-r, r])
            mot_rod.set_offsets((ps[frame][0], ps[frame][1]))
            if frame >= 3:
                velocity.set_data(vts[:frame], vps[:frame])
            if frame >= 4:
                acceleration.set_data(ats[:frame], aps[:frame])
            displacement.set_data(ts[:frame], length_with_yolk[:frame])

        self.ani = FuncAnimation(fig, update, frames=np.linspace(1, 2499, 2499), interval=50)
        self.anit = FuncAnimation(figt, update, frames=np.linspace(1, 2499, 2499), interval=50)
        self.canvas.draw()
        self.canvast.draw()


if __name__ == "__main__":
    root = tk.Tk()
    root.title('Piston Simulation')
    MainWindow(root).grid(row=0, column=0)
    root.mainloop()
