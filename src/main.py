import argparse
import sys

import numpy as np

# custom imports
from graphics import StationaryEigFigure
from hamiltonians import SingleParticleHamiltonian


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--use_density", action='store_true', default=False,
                        help="Display the probability density of the stationary states instead.")
    parser.add_argument("--save_gif", action='store_true', default=False,
                        help="Save the animation loop as a GIF.")
    parser.add_argument("--dark_mode", action='store_true', default=False,
                        help="Toggle dark mode for plot")
    args = parser.parse_args()

    # for translating so that special functions can be parsed into corresponding numpy functions
    func_dict = {
        "sin": "np.sin",
        "cos": "np.cos",
        "sinh": "np.sinh",
        "cosh": "np.cosh",
        "tanh": "np.tanh",
        "tan": "np.tan",
        "log": "np.log",
        "log10": "np.log10",
        "ln": "np.log",
        "exp": "np.exp",
        "e": "(np.e)"
    }

    # get the user's input directly instead of through argparse and also do argument checking
    user_entered_potential = input("Enter the potential (as a function of x): ")
    
    # replace with corresponding numpy functions
    for keyword in func_dict:
        user_entered_potential = user_entered_potential.replace(keyword, func_dict[keyword])

    # now get the left and right endpoints of the domain of solution
    while True:
        try:
            left_end = input("Enter the left endpoint of the interval: ")
            left_end = float(left_end)
            break
        except ValueError:
            print("Error: left endpoint could not be converted to a real number, please try again.")
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again.")
    while True:
        try:
            right_end = input("Enter the right endpoint of the interval: ")
            right_end = float(right_end)
            break
        except ValueError:
            print("Error: right endpoint could not be converted to a real number, please try again.")
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again.")

    domain_interval = [left_end, right_end]

    # need to verify that the potential works over the domain with no singularities
    pot_func = lambda x: eval(user_entered_potential)
    try:
        pot_array = pot_func(np.linspace(*domain_interval, 100))
    except NameError:
        raise NameError("Please check the spelling/formatting in the entered potential. Exiting..")
    except Exception as e:
        print("An error occurred:", str(e), ". Please double check the spelling/formatting in the entered potential before trying again. Exiting..")
        sys.exit(1)

    # input checking for number of eigenstates to find, as well as number of grid points in discretization
    while True:
        try:
            num_eigenstates_to_find = input("Enter the amount of eigenstates you want to solve for: ")
            num_grdpts = input("Enter the number of grid points in discretization: ")
            num_eigenstates_to_find = int(num_eigenstates_to_find)
            num_grdpts = int(num_grdpts)
            if not isinstance(num_eigenstates_to_find, int):
                raise ValueError("Provided number of eigenstates to find is not an integer.")
            elif not isinstance(num_grdpts, int):
                raise ValueError("Provided number of discretization (grid-) points is not an integer.")
            elif num_eigenstates_to_find > num_grdpts:
                raise ValueError(f"The maximum number of eigenstates that can be found is the number of gridpoints: {num_grdpts}.")
            break
        except ValueError as ve:
            print(ve)
            sys.exit(1)
        except Exception as e:
            print("An error occurred:", str(e), ". Please try again. Exiting..")
            sys.exit(1)

    announcement = f"""
    Solving for {num_eigenstates_to_find} eigenstates for the potential: {user_entered_potential} on {domain_interval} with {num_grdpts} grid points..
    """
    print(announcement)
    
    model = SingleParticleHamiltonian(V=pot_func, domain=domain_interval, num_points=num_grdpts)
    es, eigs = model.tise_solve(num_eigenstates_to_find)

    if args.use_density:
        print("Option '--use_density' selected. Displaying eigenstate probability densities (Born rule).")

    fig = StationaryEigFigure(energies=es, eigenstates=eigs, Vname=user_entered_potential, 
                              domain=domain_interval, density=args.use_density,
                              num_points=num_grdpts, darkmode=args.dark_mode)

    fig.begin_interactive_animate()

    # give the option to user whether to save the animation or not
    if args.save_gif:
        print("Option '--save_gif' selected. Saving the animation as a GIF.")
        fig.make_gif()
    else:
        while True:
            save_gif_input = input("Do you want to save the animation as a GIF? Answer (y/n): ")
            if 'n' in save_gif_input:
                print("Exiting without making GIF..")
                break
            elif 'y' in save_gif_input:
                print("Affirmative. Making GIF..")
                fig.make_gif()
                break
            else:
                print("Invalid option given. Please try again.")