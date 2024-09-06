import math
import sympy
from tkinter import *
from tkinter import messagebox
from tkinter import ttk

BACKGROUND = "#DBE2EF"
FONT_NAME = "Arial"
FONT_COLOR = "#393E46"


class Normal_shock:
    def __init__(self, m1, gamma):
        self.m1 = m1
        self.gamma = gamma

    def down_stream_mach(self):
        m2 = math.sqrt((1 + ((self.gamma - 1) / 2) * self.m1 ** 2) / (
                self.gamma * self.m1 ** 2 - ((self.gamma - 1) / 2)))
        return m2

    def pressure_ratio(self):
        p2_p1 = 1 + ((2 * self.gamma) / (self.gamma + 1)) * (self.m1 ** 2 - 1)
        return p2_p1

    def density_ratio(self):
        r2_r1 = ((self.gamma + 1) * self.m1 ** 2) / (2 + (self.gamma - 1) * self.m1 ** 2)
        return r2_r1

    def temperature_ratio(self):
        t2_t1 = self.pressure_ratio() / self.density_ratio()
        return t2_t1

    @staticmethod
    def total_pressure_ratio(gamma, mach):
        p0_p = (1 + ((gamma - 1) / 2) * mach ** 2) ** (gamma / (gamma - 1))
        return p0_p

    @staticmethod
    def total_temperature_ratio(gamma, mach):
        t0_t = (1 + ((gamma - 1) / 2) * mach ** 2)
        return t0_t

    def down_total_up_total_pressure(self):
        p02_p01 = (self.total_pressure_ratio(self.gamma, self.down_stream_mach()) / self.total_pressure_ratio(
            self.gamma, self.m1)) * self.pressure_ratio()
        return p02_p01

    def down_total_up_static_pressure(self):
        p02_p1 = self.total_pressure_ratio(self.gamma, self.down_stream_mach()) * self.pressure_ratio()
        return p02_p1


class Tools:
    @staticmethod
    def mach_pitot_calculate(gamma, p0, p):
        def mach_cal(x):
            m = sympy.symbols('m')
            func = (p0 / p) - (
                    ((((gamma + 1) ** 2 * m ** 2) / (4 * gamma * m ** 2 - 2 * (gamma - 1))) ** (
                            gamma / (gamma - 1))) * (
                            (1 - gamma + 2 * gamma * m ** 2) / (gamma + 1)))
            func_derivative = sympy.diff(func, m)
            iteration = 1
            while abs(func.subs(m, x)) >= 0.00001:
                if iteration >= 20:
                    return
                else:
                    x = x - (func.subs(m, x) / func_derivative.subs(m, x))
                    iteration += 1
            return x

        if gamma <= 0 or p0 <= 0 or p <= 0:
            messagebox.showinfo(title="Math Error", message="Please recheck the values.")
            raise ValueError

        try:
            mach_1 = math.sqrt((2 / (gamma - 1)) * ((p0 / p) ** ((gamma - 1) / gamma) - 1))
        except (ValueError, ZeroDivisionError) as e:
            print('Math Error, please recheck the values !!!')
            messagebox.showinfo(title="Math Error", message="Please recheck the values.")
            raise ValueError

        if mach_1 < 0.3:
            messagebox.showinfo(title="Incompressible Flow", message="The flow is incompressible.")
        elif mach_1 >= 1:
            for i in range(1, 6):
                mach_1 = mach_cal(i)
                if mach_1 and mach_1 > 0:
                    break
        if not mach_1:
            messagebox.showinfo(title="No Solution", message="Could not find any root after 20 iterations.")
        else:
            return mach_1

    @staticmethod
    def oblique_calculate(m1, theta, gamma):
        x = sympy.symbols('x')

        if m1 <= 0 or theta <= 0 or gamma <= 0:
            messagebox.showinfo(title="Math Error", message="Please recheck the values.")
            raise ValueError
        elif m1 < 1:
            messagebox.showinfo(title='Subsonic Flow', message="No shockwave for subsonic flow.")
            raise Exception('Subsonic Flow')

        a = (1 + ((gamma - 1) / 2) * m1 ** 2) * sympy.tan(theta)
        b = m1 ** 2 - 1
        c = (1 + ((gamma + 1) / 2) * m1 ** 2) * sympy.tan(theta)
        func = (a * x ** 3 - b * x ** 2 + c * x + 1)
        func_deri = (3 * a * x ** 2 - 2 * b * x + c)

        def calculate_all_root(tan_beta):
            iteration = 1
            while abs(func.subs(x, tan_beta)) >= 0.0001:
                if iteration >= 20:
                    break
                else:
                    tan_beta = tan_beta - func.subs(x, tan_beta) / func_deri.subs(x, tan_beta)
                    iteration += 1
            return tan_beta

        root = {calculate_all_root(sympy.tan(angle * (math.pi / 180))) for angle in range(0, 90, 10) if
                calculate_all_root(sympy.tan(angle * (math.pi / 180))) > 0}
        if len(root) == 0:
            messagebox.showinfo(title="No Solution", message="Could not find any root after 20 iterations.")
            return
        beta = sympy.atan(min(root))
        normal_mach_component = Normal_shock(m1 * math.sin(beta), gamma)
        m2 = normal_mach_component.down_stream_mach() / math.sin(beta - theta)
        p2p1 = normal_mach_component.pressure_ratio()
        r2r1 = normal_mach_component.density_ratio()
        t2t1 = normal_mach_component.temperature_ratio()
        p02p01 = (normal_mach_component.total_pressure_ratio(gamma, m2) / normal_mach_component.total_pressure_ratio(
            gamma, m1)) * normal_mach_component.pressure_ratio()
        p02p1 = normal_mach_component.total_pressure_ratio(gamma, m2) * normal_mach_component.pressure_ratio()
        return {"m2": m2, "beta": beta, "p2p1": p2p1, "r2r1": r2r1, "t2t1": t2t1, "p02p01": p02p01, "p02p1": p02p1}

    @staticmethod
    def expand_calculate(m1, theta, gamma):
        if m1 <= 0 or theta <= 0 or gamma <= 0:
            messagebox.showinfo(title="Math Error", message="Please recheck the values.")
            raise ValueError
        elif m1 < 1:
            messagebox.showinfo(title='Subsonic Flow', message="No shockwave for subsonic flow.")
            raise Exception('Subsonic Flow')

        prandtl_var = sympy.symbols('mach_var')
        m2 = sympy.symbols('m2')
        nuy = sympy.sqrt((gamma + 1) / (gamma - 1)) * sympy.atan(
            sympy.sqrt(((gamma - 1) / (gamma + 1)) * (prandtl_var ** 2 - 1))) - sympy.atan(
            sympy.sqrt(prandtl_var ** 2 - 1))
        func = theta + nuy.subs(prandtl_var, m1) - nuy.subs(prandtl_var, m2)
        func_deri = sympy.diff(func, m2)
        mach_solution = m1 + 0.5
        iteration = 1
        while abs(func.subs(m2, mach_solution)) > 0.0001:
            if iteration >= 20:
                messagebox.showinfo(title="No Solution", message="Could not find any root after 20 iterations.")
                return
            else:
                mach_solution = mach_solution - func.subs(m2, mach_solution) / func_deri.subs(m2, mach_solution)
                iteration += 1

        p2_p1 = Normal_shock.total_pressure_ratio(gamma, m1) / Normal_shock.total_pressure_ratio(gamma, mach_solution)
        t2_t1 = Normal_shock.total_temperature_ratio(gamma, m1) / Normal_shock.total_temperature_ratio(gamma,
                                                                                                       mach_solution)
        return {"m2": mach_solution, "p2p1": p2_p1, "t2t1": t2_t1}

    @staticmethod
    def supersonic_airfoil_2D(m1, alpha, delta, gamma, thickness, chord):
        class surface_detail:
            def __init__(self):
                self.name = None
                self.shock_type = None

        m1 = float(m1)
        alpha = float(alpha)
        gamma = float(gamma)
        delta = float(delta)
        thickness = float(thickness)
        chord = float(chord)
        description = {"top_front": surface_detail(), "top_back": surface_detail(), "bottom_front": surface_detail(),
                       "bottom_back": surface_detail()}
        if abs(alpha) >= 90:
            messagebox.showinfo('Error', 'Attack angle must be in the range of (-90;90)')
            raise ValueError
        elif abs(delta) >= 90:
            messagebox.showinfo('Error', 'Ramp angle must be in the range of (-90;90)')
            raise ValueError
        else:
            if alpha >= 0:
                description["top_front"].name = "3"
                description["top_back"].name = "5"
                description["bottom_front"].name = "2"
                description["bottom_back"].name = "4"
            else:
                description["top_front"].name = "2"
                description["top_back"].name = "4"
                description["bottom_front"].name = "3"
                description["bottom_back"].name = "5"

        # calculate
        # top front
        alpha = abs(alpha)
        delta = abs(delta)
        if alpha < delta:
            # top front oblique
            top_front_surface = Tools.oblique_calculate(gamma=gamma, m1=m1, theta=delta - alpha)
            top_front_surface_p_p1 = top_front_surface["p2p1"]
            top_front_surface_mach = top_front_surface["m2"]
            description["top_front"].shock_type = "oblique shock"
        elif alpha == delta:
            # top front no shockwave
            top_front_surface_p_p1 = 1
            top_front_surface_mach = m1
            description["top_front"].shock_type = "no shockwave"
        else:
            # top front expansion
            top_front_surface = Tools.expand_calculate(gamma=gamma, m1=m1, theta=alpha - delta)
            top_front_surface_p_p1 = top_front_surface["p2p1"]
            top_front_surface_mach = top_front_surface["m2"]
            description["top_front"].shock_type = "expansion wave"
        # top back
        top_back_surface = Tools.expand_calculate(gamma=gamma, m1=top_front_surface_mach, theta=2 * delta)
        top_back_surface_p_p1 = top_back_surface["p2p1"] * top_front_surface_p_p1
        top_back_surface_mach = top_back_surface["m2"]
        description["top_back"].shock_type = "expansion wave"

        # bottom front
        bottom_front_surface = Tools.oblique_calculate(gamma=gamma, m1=m1, theta=delta + alpha)
        bottom_front_surface_p_p1 = bottom_front_surface["p2p1"]
        bottom_front_surface_mach = bottom_front_surface["m2"]
        description["bottom_front"].shock_type = "oblique shock"

        # bottom back
        bottom_back_surface = Tools.expand_calculate(gamma=gamma, m1=bottom_front_surface_mach, theta=2 * delta)
        bottom_back_surface_p_p1 = bottom_back_surface["p2p1"] * bottom_front_surface_p_p1
        bottom_back_surface_mach = bottom_back_surface["m2"]
        description["bottom_back"].shock_type = "expansion wave"

        # drag lift coefficient
        cl = (1 / (2 * math.cos(delta) * 0.5 * gamma * m1 ** 2)) * (
                (bottom_front_surface_p_p1 - top_back_surface_p_p1) * math.cos(delta + alpha) + (
                bottom_back_surface_p_p1 - top_front_surface_p_p1) * math.cos(delta - alpha))
        cd = (1 / (2 * math.cos(delta) * 0.5 * gamma * m1 ** 2)) * (
                (bottom_front_surface_p_p1 - top_back_surface_p_p1) * math.sin(delta + alpha) + (
                -bottom_back_surface_p_p1 + top_front_surface_p_p1) * math.sin(delta - alpha))

        return {"cl": cl, "cd": cd, "description": description}

    @staticmethod
    def supersonic_airfoil_3D(m1, alpha, delta, gamma, thickness, chord, lamda):
        m1_equi = m1 * math.sqrt(1 - (math.sin(lamda) * math.cos(alpha)) ** 2)
        alpha_equi = math.atan(math.tan(alpha) / math.cos(lamda))

        [cl2D, cd2D, description] = Tools.supersonic_airfoil_2D(m1_equi, alpha_equi, delta, gamma, thickness,
                                                                chord).values()
        cl3D = float(cl2D) * (1 - (math.sin(lamda) * math.cos(alpha)) ** 2)
        cd3D = float(cd2D) * math.cos(lamda) * (1 - (math.sin(lamda) * math.cos(alpha)) ** 2)

        return {"cl": cl3D, "cd": cd3D, "description": description}


class Gui:
    # Gui creation
    def __init__(self):
        window = Tk()
        window.title("Compressible Flow Calculator")
        window.config(bg=BACKGROUND)
        window.minsize(670, 600)
        window.geometry("670x600")

        style = ttk.Style(window)
        style.theme_create("DLYN", parent="alt", settings={
            ".": {
                "configure": {
                    "background": BACKGROUND,  # All colors except for active tab-button
                    "font": 'red'
                }
            },
            "TNotebook.Tab": {
                "configure": {
                    "background": BACKGROUND,  # Color of non-selected tab-button
                    "padding": [5, 2],
                    # [space between text and horizontal tab-button border, space between text and vertical tab_button border]
                    "font": "white"
                },
                "map": {
                    "background": [("selected", '#FAF7F0')],  # Color of active tab
                    "expand": [("selected", [1, 1, 1, 0])]  # [expanse of text]
                }
            }

        })
        style.theme_use("DLYN")
        style.configure('lefttab.TNotebook', tabposition='wn')

        # Window creation
        tabControl = ttk.Notebook(window,style='lefttab.TNotebook')
        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)
        tab3 = ttk.Frame(tabControl)
        tab4 = ttk.Frame(tabControl)
        tab5 = ttk.Frame(tabControl)

        tabControl.add(tab1, text=f'{"Pitot Mach": ^25s}')
        tabControl.add(tab2, text=f'{"Normal Shock": ^20s}')
        tabControl.add(tab3, text=f'{"Oblique Shock": ^20s}')
        tabControl.add(tab4, text=f'{"Expansion Fan": ^20s}')
        tabControl.add(tab5, text=f'{"Supersonic Airfoil": ^20s}')

        tabControl.pack(side=TOP,fill=BOTH,expand=True)

        # Tab1 creation
        canvas1 = Canvas(tab1, width=460, height=200, bg=BACKGROUND, highlightthickness=0)
        pitot_img = PhotoImage(file="pitot-tube.png")
        canvas1.create_image(230, 100, image=pitot_img)
        canvas1.create_text(90, 100, text="Po", fill=FONT_COLOR, font=(FONT_NAME, 15, "bold"))
        canvas1.create_text(140, 100, text="P", fill=FONT_COLOR, font=(FONT_NAME, 15, "bold"))
        canvas1.pack()

        down_frame_mach = Frame(tab1, bg=BACKGROUND)
        down_frame_mach.pack()
        Label(down_frame_mach, text="γ", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1, column=0,
                                                                                                  padx=10)
        self.gamma_text = Entry(down_frame_mach, width=20)
        self.gamma_text.grid(row=1, column=1)
        Label(down_frame_mach, text="Po (Pa)", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND).grid(row=2, column=0,
                                                                                                   padx=10)
        self.P0_text = Entry(down_frame_mach, width=20)
        self.P0_text.grid(row=2, column=1)
        Label(down_frame_mach, text="P (Pa)", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND).grid(row=3, column=0,
                                                                                                  padx=10)
        self.P_text = Entry(down_frame_mach, width=20)
        self.P_text.grid(row=3, column=1)
        mach_calculation_button = Button(down_frame_mach, text="Calculate", command=self.tab1_button,
                                         highlightthickness=0,
                                         font=(FONT_NAME, 12,))
        mach_calculation_button.grid(row=4, column=0, columnspan=2, pady=10)
        self.mach_label = Label(down_frame_mach, text="M = ", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.mach_label.grid(row=5, column=0, columnspan=2)

        # Tab2 creation
        up_frame_normal = Frame(tab2, bg=BACKGROUND)
        up_frame_normal.pack(pady=20)
        down_frame_normal = Frame(tab2, bg=BACKGROUND)
        down_frame_normal.pack()

        Label(up_frame_normal, text="γ", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND).grid(row=0, column=0,
                                                                                                  padx=10)
        self.gamma_normal_text = Entry(up_frame_normal, width=20)
        self.gamma_normal_text.grid(row=0, column=1)
        Label(up_frame_normal, text="M1", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1, column=0,
                                                                                                   padx=10)
        self.M1_text_normal = Entry(up_frame_normal, width=20)
        self.M1_text_normal.grid(row=1, column=1)
        normal_shock_button = Button(up_frame_normal, text="Calculate", command=self.tab2_button, highlightthickness=0,
                                     font=(FONT_NAME, 12,))
        normal_shock_button.grid(row=2, column=0, columnspan=2, pady=10)

        Label(down_frame_normal, text="M2", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=0, column=0, )
        Label(down_frame_normal, text="p2/p1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1,
                                                                                                         column=0, )
        Label(down_frame_normal, text="ρ2/ρ1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=2,
                                                                                                         column=0, )
        Label(down_frame_normal, text="T2/T1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=3,
                                                                                                         column=0, )
        Label(down_frame_normal, text="p02/p01", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=4,
                                                                                                           column=0)
        Label(down_frame_normal, text="p02/p1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=5,
                                                                                                          column=0)
        for row in range(0, len(down_frame_normal.winfo_children())):
            Label(down_frame_normal, text=" = ", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=row,
                                                                                                           column=1)
        self.M2_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.M2_label_normal.grid(row=0, column=2, )
        self.p2_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.p2_label_normal.grid(row=1, column=2, )
        self.r2_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.r2_label_normal.grid(row=2, column=2, )
        self.T2_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.T2_label_normal.grid(row=3, column=2, )
        self.p02_p01_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                          bg=BACKGROUND)
        self.p02_p01_label_normal.grid(row=4, column=2)
        self.p02_p1_label_normal = Label(down_frame_normal, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                         bg=BACKGROUND)
        self.p02_p1_label_normal.grid(row=5, column=2)

        # Tab3 creation
        canvas3 = Canvas(tab3, width=460, height=200, bg=BACKGROUND, highlightthickness=0)
        oblique_img = PhotoImage(file="rsz_oblique.png")
        canvas3.create_image(230, 100, image=oblique_img)
        canvas3.create_text(230, 20, text="Week shock calculation", fill=FONT_COLOR, font=(FONT_NAME, 12, "bold"))
        canvas3.pack()

        up_frame_oblique = Frame(tab3, bg=BACKGROUND)
        up_frame_oblique.pack(pady=10)
        down_frame_oblique = Frame(tab3, bg=BACKGROUND)
        down_frame_oblique.pack(pady=5)

        gamma_label_oblique_label = Label(up_frame_oblique, text="γ", font=(FONT_NAME, 15), fg=FONT_COLOR,
                                          bg=BACKGROUND)
        gamma_label_oblique_label.grid(row=1, column=0, padx=10)
        self.gamma_oblique_text = Entry(up_frame_oblique, width=20)
        self.gamma_oblique_text.grid(row=1, column=1)
        theta_label = Label(up_frame_oblique, text="θ (degree)", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND)
        theta_label.grid(row=2, column=0, padx=10)
        self.theta_text = Entry(up_frame_oblique, width=20)
        self.theta_text.grid(row=2, column=1)
        M1_label_oblique = Label(up_frame_oblique, text="M1", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND)
        M1_label_oblique.grid(row=3, column=0, padx=10)
        self.M1_text_oblique = Entry(up_frame_oblique, width=20)
        self.M1_text_oblique.grid(row=3, column=1)
        beta_calculation_button = Button(up_frame_oblique, text="Calculate", command=self.tab3_button,
                                         highlightthickness=0,
                                         font=(FONT_NAME, 12,))
        beta_calculation_button.grid(row=5, column=0, columnspan=2, pady=10)

        Label(down_frame_oblique, text="β", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=0, column=0, )
        Label(down_frame_oblique, text="M2", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1,
                                                                                                       column=0, )
        Label(down_frame_oblique, text="p2/p1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=2,
                                                                                                          column=0, )
        Label(down_frame_oblique, text="ρ2/ρ1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=3,
                                                                                                          column=0, )
        Label(down_frame_oblique, text="T2/T1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=4,
                                                                                                          column=0, )
        Label(down_frame_oblique, text="p02/p01", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=5,
                                                                                                            column=0)
        Label(down_frame_oblique, text="p02/p1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=6,
                                                                                                           column=0)
        for row in range(0, len(down_frame_oblique.winfo_children())):
            Label(down_frame_oblique, text=" = ", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=row,
                                                                                                            column=1)

        self.beta_label = Label(down_frame_oblique, text="", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND)
        self.beta_label.grid(row=0, column=2)
        self.M2_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.M2_label_oblique.grid(row=1, column=2, )
        self.p2_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.p2_label_oblique.grid(row=2, column=2, )
        self.r2_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.r2_label_oblique.grid(row=3, column=2, )
        self.T2_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND)
        self.T2_label_oblique.grid(row=4, column=2, )
        self.p02_p01_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                           bg=BACKGROUND)
        self.p02_p01_label_oblique.grid(row=5, column=2)
        self.p02_p1_label_oblique = Label(down_frame_oblique, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                          bg=BACKGROUND)
        self.p02_p1_label_oblique.grid(row=6, column=2)

        # tab4 creation
        canvas4 = Canvas(tab4, width=460, height=250, bg=BACKGROUND, highlightthickness=0)
        expansion_img = PhotoImage(file="expansion-fan.png")
        canvas4.create_image(230, 125, image=expansion_img)
        canvas4.pack()

        up_frame_expansion = Frame(tab4, bg=BACKGROUND)
        up_frame_expansion.pack(pady=10)
        down_frame_expansion = Frame(tab4, bg=BACKGROUND)
        down_frame_expansion.pack()

        gamma_label_expansion_label = Label(up_frame_expansion, text="γ", font=(FONT_NAME, 15), fg=FONT_COLOR,
                                            bg=BACKGROUND)
        gamma_label_expansion_label.grid(row=0, column=0, padx=10)
        self.gamma_expansion_text = Entry(up_frame_expansion, width=20)
        self.gamma_expansion_text.grid(row=0, column=1)
        theta_expansion_label = Label(up_frame_expansion, text="θ (degree)", font=(FONT_NAME, 15), fg=FONT_COLOR,
                                      bg=BACKGROUND)
        theta_expansion_label.grid(row=1, column=0, padx=10)
        self.theta_expansion_text = Entry(up_frame_expansion, width=20)
        self.theta_expansion_text.grid(row=1, column=1)
        M1_label_expansion = Label(up_frame_expansion, text="M1", font=(FONT_NAME, 15), fg=FONT_COLOR, bg=BACKGROUND)
        M1_label_expansion.grid(row=2, column=0, padx=10)
        self.M1_text_expansion = Entry(up_frame_expansion, width=20)
        self.M1_text_expansion.grid(row=2, column=1)
        m2_calculation_button = Button(up_frame_expansion, text="Calculate", command=self.tab4_button,
                                       highlightthickness=0,
                                       font=(FONT_NAME, 12,))
        m2_calculation_button.grid(row=3, column=0, columnspan=2, pady=10)

        Label(down_frame_expansion, text="M2", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=0,
                                                                                                         column=0, )
        Label(down_frame_expansion, text="p2/p1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1,
                                                                                                            column=0, )
        Label(down_frame_expansion, text="T2/T1", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=2,
                                                                                                            column=0, )
        for row in range(0, len(down_frame_expansion.winfo_children())):
            Label(down_frame_expansion, text=" = ", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=row,
                                                                                                              column=1)

        self.M2_label_expansion = Label(down_frame_expansion, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                        bg=BACKGROUND)
        self.M2_label_expansion.grid(row=0, column=2, )
        self.p2_label_expansion = Label(down_frame_expansion, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                        bg=BACKGROUND)
        self.p2_label_expansion.grid(row=1, column=2, )
        self.T2_label_expansion = Label(down_frame_expansion, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                        bg=BACKGROUND)
        self.T2_label_expansion.grid(row=2, column=2, )

        # tab5 creation
        canvas5 = Canvas(tab5, width=460, height=250, bg=BACKGROUND, highlightthickness=0)
        double_wedge_img = PhotoImage(file="double-wedge.png")
        canvas5.create_image(230, 125, image=double_wedge_img)
        canvas5.pack()

        up_frame_supersonic_airfoil = Frame(tab5, bg=BACKGROUND)
        up_frame_supersonic_airfoil.pack(pady=5)
        down_frame_supersonic_airfoil = Frame(tab5, bg=BACKGROUND)
        down_frame_supersonic_airfoil.pack()

        gamma_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="γ", font=(FONT_NAME, 15),
                                               fg=FONT_COLOR,
                                               bg=BACKGROUND)
        gamma_label_supersonic_airfoil.grid(row=0, column=0, padx=10)
        self.gamma_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.gamma_supersonic_airfoil.grid(row=0, column=1)
        M1_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="M1", font=(FONT_NAME, 15), fg=FONT_COLOR,
                                            bg=BACKGROUND)
        M1_label_supersonic_airfoil.grid(row=1, column=0, padx=10)
        self.M1_text_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.M1_text_supersonic_airfoil.grid(row=1, column=1)

        thickness_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="Thickness", font=(FONT_NAME, 15),
                                                   fg=FONT_COLOR, bg=BACKGROUND)
        thickness_label_supersonic_airfoil.grid(row=2, column=0, padx=10)
        self.thickness_text_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.thickness_text_supersonic_airfoil.grid(row=2, column=1)

        chord_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="Chord length", font=(FONT_NAME, 15),
                                               fg=FONT_COLOR, bg=BACKGROUND)
        chord_label_supersonic_airfoil.grid(row=3, column=0, padx=10)
        self.chord_text_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.chord_text_supersonic_airfoil.grid(row=3, column=1)

        delta_supersonic_airfoil_label = Label(up_frame_supersonic_airfoil, text="δ (degree)", font=(FONT_NAME, 15),
                                               fg=FONT_COLOR,
                                               bg=BACKGROUND)
        delta_supersonic_airfoil_label.grid(row=0, column=3, padx=10)
        self.delta_supersonic_airfoil_text = Entry(up_frame_supersonic_airfoil, width=17)
        self.delta_supersonic_airfoil_text.grid(row=0, column=4)

        alpha_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="α (degree)", font=(FONT_NAME, 15),
                                               fg=FONT_COLOR, bg=BACKGROUND)
        alpha_label_supersonic_airfoil.grid(row=1, column=3, padx=10)
        self.alpha_text_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.alpha_text_supersonic_airfoil.grid(row=1, column=4)
        lambda_label_supersonic_airfoil = Label(up_frame_supersonic_airfoil, text="Λ (degree)", font=(FONT_NAME, 15),
                                                fg=FONT_COLOR, bg=BACKGROUND)
        lambda_label_supersonic_airfoil.grid(row=2, column=3, padx=10)
        self.lambda_text_supersonic_airfoil = Entry(up_frame_supersonic_airfoil, width=17)
        self.lambda_text_supersonic_airfoil.grid(row=2, column=4)
        supersonic_airfoil_calculation_button = Button(up_frame_supersonic_airfoil, text="Calculate",
                                                       command=self.tab5_button,
                                                       highlightthickness=0,
                                                       font=(FONT_NAME, 12,))
        supersonic_airfoil_calculation_button.grid(row=4, column=0, columnspan=4, pady=10)

        Label(down_frame_supersonic_airfoil, text="Cl", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=0,
                                                                                                                  column=0, )
        Label(down_frame_supersonic_airfoil, text="Cd", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(row=1,
                                                                                                                  column=0, )
        for row in range(0, len(down_frame_supersonic_airfoil.winfo_children())):
            Label(down_frame_supersonic_airfoil, text=" = ", font=(FONT_NAME, 15,), fg=FONT_COLOR, bg=BACKGROUND).grid(
                row=row,
                column=1)

        self.cl_supersonic_airfoil = Label(down_frame_supersonic_airfoil, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                           bg=BACKGROUND)
        self.cl_supersonic_airfoil.grid(row=0, column=2)
        self.cd_supersonic_airfoil = Label(down_frame_supersonic_airfoil, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                           bg=BACKGROUND)
        self.cd_supersonic_airfoil.grid(row=1, column=2)
        self.shock_wave_detail = Label(down_frame_supersonic_airfoil, text="", font=(FONT_NAME, 15,), fg=FONT_COLOR,
                                       bg=BACKGROUND, justify=LEFT)
        self.shock_wave_detail.grid(row=2, column=0, columnspan=3)
        window.mainloop()

    def tab1_button(self):
        gamma = float(self.gamma_text.get())
        p0 = float(self.P0_text.get())
        p = float(self.P_text.get())
        self.mach_label.config(text=("M = {:.4f}".format(Tools.mach_pitot_calculate(gamma, p0, p))))

    def tab2_button(self):
        m1 = self.M1_text_normal.get()
        gamma = self.gamma_normal_text.get()
        if float(m1) <= 0 or float(m1) <= 0:
            messagebox.showinfo(title='Math Error', message="Please recheck the values.")
            raise ValueError
        elif float(m1) < 1:
            messagebox.showinfo(title='Subsonic Flow', message="No shockwave for subsonic flow.")
            raise Exception('Subsonic Flow')
        else:
            normal_shock = Normal_shock(m1=float(m1), gamma=float(gamma))
            self.M2_label_normal.config(text=("{:.4f}".format(normal_shock.down_stream_mach())))
            self.p2_label_normal.config(text=("{:.4f}".format(normal_shock.pressure_ratio())))
            self.r2_label_normal.config(text=("{:.4f}".format(normal_shock.density_ratio())))
            self.T2_label_normal.config(text=("{:.4f}".format(normal_shock.temperature_ratio())))
            self.p02_p01_label_normal.config(text=("{:.4f}".format(normal_shock.down_total_up_total_pressure())))
            self.p02_p1_label_normal.config(text=("{:.4f}".format(normal_shock.down_total_up_static_pressure())))

    def tab3_button(self):
        m1 = float(self.M1_text_oblique.get())
        theta = float(self.theta_text.get()) * (math.pi / 180)
        gamma = float(self.gamma_oblique_text.get())
        oblique_solver = Tools.oblique_calculate(m1, theta, gamma)
        self.beta_label.config(text=("{:.4f}°".format(oblique_solver["beta"] * (180 / math.pi))))
        self.M2_label_oblique.config(text=("{:.4f}".format(oblique_solver["m2"])))
        self.p2_label_oblique.config(text=("{:.4f}".format(oblique_solver["p2p1"])))
        self.r2_label_oblique.config(text=("{:.4f}".format(oblique_solver["r2r1"])))
        self.T2_label_oblique.config(text=("{:.4f}".format(oblique_solver["t2t1"])))
        self.p02_p01_label_oblique.config(text=("{:.4f}".format(oblique_solver["p02p01"])))
        self.p02_p1_label_oblique.config(text=("{:.4f}".format(oblique_solver["p02p1"])))

    def tab4_button(self):
        m1 = float(self.M1_text_expansion.get())
        theta = float(self.theta_expansion_text.get()) * (math.pi / 180)
        gamma = float(self.gamma_expansion_text.get())
        expansion_solve = Tools.expand_calculate(m1, theta, gamma)
        self.M2_label_expansion.config(text=("{:.4f}".format(expansion_solve["m2"])))
        self.p2_label_expansion.config(text=("{:.4f}".format(expansion_solve["p2p1"])))
        self.T2_label_expansion.config(text=("{:.4f}".format(expansion_solve["t2t1"])))

    def tab5_button(self):
        m1 = float(self.M1_text_supersonic_airfoil.get())
        delta = float(self.delta_supersonic_airfoil_text.get()) * (math.pi / 180)
        gamma = float(self.gamma_supersonic_airfoil.get())
        thickness = float(self.thickness_text_supersonic_airfoil.get())
        chord = float(self.chord_text_supersonic_airfoil.get())
        if len(self.alpha_text_supersonic_airfoil.get()) != 0:
            alpha = float(self.alpha_text_supersonic_airfoil.get()) * (math.pi / 180)
        else:
            alpha = 0
        if len(self.lambda_text_supersonic_airfoil.get()) != 0:
            lamda = float(self.lambda_text_supersonic_airfoil.get()) * (math.pi / 180)
        else:
            lamda = 0

        if lamda != 0:
            [cl, cd, description] = Tools.supersonic_airfoil_3D(m1=m1, alpha=alpha, delta=delta, lamda=lamda,
                                                                gamma=gamma, thickness=thickness, chord=chord).values()
        else:
            [cl, cd, description] = Tools.supersonic_airfoil_2D(m1=m1, alpha=alpha, delta=delta,
                                                                gamma=gamma, thickness=thickness, chord=chord).values()
        self.cl_supersonic_airfoil.config(text=("{:.4f}".format(cl)))
        self.cd_supersonic_airfoil.config(text=("{:.4f}".format(cd)))
        text = ""
        for surface in description:
            surface_name = description[surface].name
            wave_type = description[surface].shock_type
            text += f"{surface_name}   :   {wave_type}\n"
        self.shock_wave_detail.config(text=text)


if __name__ == "__main__":
    gui = Gui()
# Tools.supersonic_airfoil_3D(2, 5 * (math.pi / 180), 5 * (math.pi / 180), 1.4, 0.08749, 1, 30 * (math.pi / 180))
