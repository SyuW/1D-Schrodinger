from scipy.integrate import odeint

import argparse
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt


# set mass and hbar to 1
def solve_time_independent_schrodinger(energy_guess):
    return


def guess_energy():
    return


def create_animation():
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enter a one-dimensional potential")
    parser.add_argument("--potential", metavar="x**2", help="potential")
    parser.add_argument("-v", type=bool, metavar="N", help="verbose")

    potential = lambda x: eval(args.v)

    args = parser.parse_args()
