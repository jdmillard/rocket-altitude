'''
The entire purpose of these functions is to find the least common multple
of n numbers. lcm_mult is the one to import

 >>> NOT CURRENTLY USED <<<
'''

# define greatest common divisor function (2 numbers)
def gcd(x, y):
    '''This function implements the Euclidian algorithm
    to find G.C.D. of two numbers'''
    while(y):
        x, y = y, x % y
    return x

# define greatest least common multiple function (2 numbers)
def lcm(x,y):
    '''This function implements the Euclidian algorithm
    to find L.C.M. of two numbers'''
    out = (x*y) // gcd(x,y) # // performs integer division
    return out

# define the greatest common multiple function (n numbers)
def lcm_mult(list_in):
    '''This function uses the fact that lcm(a,b,c) = lcm(a,lcm(b,c))
    to find L.C.M. of n numbers'''
    lcm_cur = 1 # initialize the current lcm value
    for i in list_in:
        lcm_cur = lcm(lcm_cur, i)
    return lcm_cur
