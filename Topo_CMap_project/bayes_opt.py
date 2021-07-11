from skopt import gp_minimize
from validation import validate

res = gp_minimize(validate, [(1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (1.0, 10.0), (0.0, 1.0)],
                  n_initial_points=15, n_calls=30, verbose=True)
print(str(res))
with open("results.txt", "w") as file:
    file.write(res)






