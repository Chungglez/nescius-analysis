from math import sqrt, cos,sin,tan,radians, log10
import numpy as np

def sind(x):
    return sin(radians(x))
def cosd(x):
    return cos(radians(x))
def tand(x):
    return tan(radians(x))


# Safe interpolation wrapper for np function
def interp(x: float, x_list, y_list) -> float:
    # Validate inputs
    if len(x_list) != len(y_list):
        raise ValueError("x_list and y_list must have the same length.")
    if len(x_list) < 2:
        raise ValueError("x_list and y_list must contain at least two elements.")
    if not np.issubdtype(type(x), np.number):
        raise TypeError("x must be a numeric value.")
    if not np.all(np.isfinite(x_list)) or not np.all(np.isfinite(y_list)):
        raise ValueError("x_list and y_list must contain finite numbers.")

    # Check monotonicity
    if all(x_list[i] <= x_list[i + 1] for i in range(len(x_list) - 1)):
        # x_list is positively increasing
        return np.interp(x, x_list, y_list)
    elif all(x_list[i] >= x_list[i + 1] for i in range(len(x_list) - 1)):
        # x_list is negatively decreasing
        return np.interp(x, x_list[::-1], y_list[::-1])
    else:
        raise ValueError("x_list must be monotonically increasing or decreasing.")


