from decimal import Decimal, getcontext
from math import log10
import matplotlib.pyplot as plt

getcontext().prec = 500

x = input("Enter the value of beta: ")
y = input("Enter the initial value of the recurrence: ")
beta = Decimal(x)
z = Decimal(y)

tolerance = Decimal("1e-8")
q_max = int(input("Enter the maximum value of q: "))
max_iterations = 1000

q_values = []
approx_errors = []
log_approx_errors = []
iterations_needed = []

for q in range(q_max + 1):

    a = beta * (Decimal(10) ** q)
    b = beta

    # Fixed point
    phi_star = (a + (a*a +4*b).sqrt())/2

    # Taylor approximation
    approximation = a + b/a

    approximation_error = abs(phi_star - approximation)

    q_values.append(q)
    approx_errors.append(float(approximation_error))

    if approximation_error != 0:
        log_approx_errors.append(float(approximation_error.log10()))
    else:
        log_approx_errors.append (None)

    Phi = z

    phi_values = []
    errors = []
    iterations = []

    i = 0

    while abs(Phi - phi_star) > tolerance and i < max_iterations:

        phi_values.append(float(Phi))
        errors.append(float(abs(Phi - phi_star)))
        iterations.append(i)

        Phi = a + b/Phi

        i += 1

    iterations_needed.append(i)

#Graph 1 - Iterations Needed x q
plt.figure(figsize = (8,5))
plt.plot(
    q_values,
    iterations_needed,
    marker = "o",
)

plt.xlabel("q")
plt.ylabel("iterations Required for Convergence")

plt.grid(True)
plt.tight_layout()
#Graph 2 - Error x q

plt.figure(figsize = (8,5))
plt.plot(
    q_values, approx_errors,marker = "o"
)
plt.yscale("log")
plt.xlabel("q")
plt.ylabel("Approximation Error")

plt.title("Taylor Approximation Error")

plt.grid(True, which = "both")
plt.tight_layout()

#Graph 3 - Slope Analysis

theory = [-3*q - log10(float(beta))
          for q in q_values]

plt.figure(figsize = (8,5))

plt.plot(
    q_values, log_approx_errors, marker = "o", label = "Numerical"
)

plt.plot(
    q_values, theory, linestyle = "--", label = "theory (slope -3")

plt.xlabel("q")
plt.ylabel("log10 error")
plt.title("Numerical vs Theoretical Error Decay")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.show()