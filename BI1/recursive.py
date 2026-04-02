def factorial(n):
    if n == 0:
        return 1  # Base case
    else:
        print(f"Value of n is: {n}")
        print("And the value of the factorial here is: ", n * factorial(n - 1))
    
print(factorial(5))