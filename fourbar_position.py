import math
import numpy as np
import matplotlib.pyplot as plt


def k1(a, d):
    return d / a


def k2(c, d):
    return d / c


def k3(a, b, c, d):
    return (a**2 - b**2 + c**2 + d**2) / (2 * a * c)


def k4(b, d):
    return d / b


def k5(a, b, c, d):
    return (c**2 - d**2 - a**2 - b**2) / (2 * a * b)


def A(a, b, c, d, theta2):
    return math.cos(theta2) - k1(a, d) - k2(c, d) * math.cos(theta2) + k3(a, b, c, d)


def B(theta2):
    return -2 * math.sin(theta2)


def C(a, b, c, d, theta2):
    return k1(a, d) - (k2(c, d) + 1) * math.cos(theta2) + k3(a, b, c, d)


def D(a, b, c, d, theta2):
    return math.cos(theta2) - k1(a, d) + k4(b, d) * math.cos(theta2) + k5(a, b, c, d)


def E(theta2):
    return -2 * math.sin(theta2)


def F(a, b, c, d, theta2):
    return k1(a, d) + (k4(b, d) - 1) * math.cos(theta2) + k5(a, b, c, d)


def theta_4(a, b, c, d, theta2, circuit="open"):
    root = B(theta2) ** 2 - 4 * A(a, b, c, d, theta2) * C(a, b, c, d, theta2)
    sign = -1 if circuit == "open" else 1
    return 2 * math.atan((-B(theta2) + sign * math.sqrt(root)) / (2 * A(a, b, c, d, theta2)))


def theta_3(a, b, c, d, theta2, circuit="open"):
    root = E(theta2) ** 2 - 4 * D(a, b, c, d, theta2) * F(a, b, c, d, theta2)
    sign = -1 if circuit == "open" else 1
    return 2 * math.atan((-E(theta2) + sign * math.sqrt(root)) / (2 * D(a, b, c, d, theta2)))


def transmission_angle(theta3, theta4):
    angle = math.fabs(theta3 - theta4)
    if angle > (math.pi / 2):
        return math.pi - angle
    return angle



def omega_3(a, b, theta2, theta3, theta4, omega2):
    return (a * omega2 / b) * math.sin(theta4 - theta2) / math.sin(theta3 - theta4)
 
 
def omega_4(a, c, theta2, theta3, theta4, omega2):
    return (a * omega2 / c) * math.sin(theta2 - theta3) / math.sin(theta4 - theta3)
 
 
def velocity_A(a, theta2, omega2):
    vx = -a * omega2 * math.sin(theta2)
    vy = a * omega2 * math.cos(theta2)
    return vx, vy
 
 
def velocity_B(c, theta4, omega4):
    vx = -c * omega4 * math.sin(theta4)
    vy = c * omega4 * math.cos(theta4)
    return vx, vy

def velocity_P(p, theta4, omega4):
    vx = -p * omega4 * math.sin(theta4)
    vy = p * omega4 * math.cos(theta4)
    return vx, vy


def rotate_velocity_vector(vx, vy, offset_deg):
    offset_rad = math.radians(offset_deg)
    cos_offset = math.cos(offset_rad)
    sin_offset = math.sin(offset_rad)

    rotated_vx = vx * cos_offset - vy * sin_offset
    rotated_vy = vx * sin_offset + vy * cos_offset
    return rotated_vx, rotated_vy
 
 
def velocity_BA(b, theta3, omega3):
    vx = -b * omega3 * math.sin(theta3)
    vy = b * omega3 * math.cos(theta3)
    return vx, vy
 
 
def speed(vx, vy):
    return math.hypot(vx, vy)
 
 
 
def slider_d_dot(a, b, theta2, theta3, omega2, omega3):
    return -a * omega2 * math.sin(theta2) + b * omega3 * math.sin(theta3)
 
 
def slider_velocity_B(d_dot):
    return d_dot, 0.0





a = 2.5
b = 3.5
c = 3.9
d = 4.638
p = 6
omega2 = 7.15

THETA2_MIN_DEG = -75.2965724
THETA2_MAX_DEG = 11.836
STEP_DEG = 1
GROUND_LINK_OFFSET_ANGLE_DEG = -14.0362435

theta2_deg = np.arange(THETA2_MIN_DEG, THETA2_MAX_DEG + STEP_DEG, STEP_DEG)


theta3_deg = np.array([math.degrees(theta_3(a, b, c, d, math.radians(t))) for t in theta2_deg])
theta3_deg = theta3_deg + GROUND_LINK_OFFSET_ANGLE_DEG
theta4_deg = np.array([math.degrees(theta_4(a, b, c, d, math.radians(t))) for t in theta2_deg])
theta4_deg = theta4_deg + GROUND_LINK_OFFSET_ANGLE_DEG
mu_deg = np.array([
    math.degrees(transmission_angle(math.radians(t3), math.radians(t4)))
    for t3, t4 in zip(theta3_deg, theta4_deg)
])
speed_A = np.array([speed(*velocity_A(a, math.radians(t), omega2)) for t in theta2_deg])
speed_P = np.array([speed(*velocity_P(p, theta_4(a, b, c, d, math.radians(t)), omega_4(a, c, math.radians(t), theta_3(a, b, c, d, math.radians(t)), theta_4(a, b, c, d, math.radians(t)), omega2))) for t in theta2_deg])
vel_Py = []
for t in theta2_deg:
    t_rad = math.radians(t)
    th3 = theta_3(a, b, c, d, t_rad)
    th4 = theta_4(a, b, c, d, t_rad)
    w4 = omega_4(a, c, t_rad, th3, th4, omega2)

    vx, vy = velocity_P(p, th4, w4)
    _, vy_rot = rotate_velocity_vector(vx, vy, GROUND_LINK_OFFSET_ANGLE_DEG)
    vel_Py.append(vy_rot)

vel_Py = np.array(vel_Py)
OFFSET_DEG = 360
theta2_deg = theta2_deg + OFFSET_DEG + GROUND_LINK_OFFSET_ANGLE_DEG

def plot_vs_theta2(y, ylabel, title=None):
    plt.figure(figsize=(9, 5))
    plt.plot(theta2_deg, y, lw=2)
    plt.xlabel(r"$\theta_2$ (deg)")
    plt.ylabel(ylabel)
    plt.xlim(theta2_deg[0], theta2_deg[-1]) 
    if title:
        plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


print(theta_4(a, b, c, d, math.radians(11.85)))

 

plot_vs_theta2(mu_deg,     "Transmission angle (deg)", "Transmission angle vs. $\\theta_2$")
plot_vs_theta2(theta3_deg, r"$\theta_3$ (deg)",        r"$\theta_3$ vs. $\theta_2$")
plot_vs_theta2(theta4_deg, r"$\theta_4$ (deg)",        r"$\theta_4$ vs. $\theta_2$")

plot_vs_theta2(speed_A, r"$|V_A|$ (in/s)", r"$|V_A|$ vs. $\theta_2$")
plot_vs_theta2(speed_P, r"$|V_C|$ (in/s)", r"$|V_C|$ vs. $\theta_2$")
plot_vs_theta2(vel_Py, r"$V_{C,y}$ (in/s)", r"$V_{C,y}$ vs. $\theta_2$")


