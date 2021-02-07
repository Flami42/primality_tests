import random
import time
import matplotlib.pyplot as plt
from varname import nameof
from random import randrange
import math
import multiprocessing

def Timer(s, flag):
    time.sleep(s)  # s seconds
    flag[0] = False  

def gcd(x,y):
    while y != 0:
        x, y = y, x%y
    return x

def is_coprime(x, y):
    return gcd(x, y) == 1


# compute (a^n)%p 
def power(a, n, p):
    # Initialize result 
    res = 1
    # Update 'a' if 'a' >= p 
    a = a % p  
    while n > 0:
        # If n is odd, multiply 'a' with result 
        if n % 2:
            res = (res * a) % p
            n = n - 1
        else:
            a = (a ** 2) % p
         # n must be even now 
            n = n // 2
             
    return res % p

    #jacobian
def legendre(a, p):
    if p < 2:
        raise ValueError('p must not be < 2')
    if (a == 0) or (a == 1):
        return a
    if a % 2 == 0:
        r = legendre(a // 2, p)
        if p * p - 1 & 8 != 0:
            r *= -1
    else:
        r = legendre(p % a, a)
        if (a - 1) * (p - 1) & 4 != 0:
            r *= -1
    return r    

######################### FACTORING #######################

#takes in input a number (variable "results" is needed just for implememtation purposes)
def factoring(n, results):
    #corner cases
    if n ==2 or n==3:
        return True
    #compute the qsquare root and test all the numbers 
    sq = int(math.sqrt(n))
    for i in range(2,sq+1):
        if n%i==0 : #if any of them (!=1) divides n, return Composite (False)
            results[n]= False
    #otherwise return Prime (True)
    results[n]= True
            
######################### FERMAT #######################

def fermat_test(n):
    # Fermat's little theorem 
    a = random.randint(2, n - 2)
    if power(a, n - 1, n) != 1:
        return False # return composite
    else:
        return True # return probably a prime

def fermat_probabilistic_test(n, k):
    # Corner cases
    if n == 1 or n == 4: 
        return False
    elif n == 2 or n == 3: # first two important primes
        return True
    #otherwise repeat fermat test k times
    else: 
        for i in range(k):
            if not fermat_test(n):
                return False
    return True

######################### MILLER RABIN #######################

def miller_rabin(n, k): 
    # Corner cases
    if n == 2 or n == 3:
        return True
    # If number is even, it's a composite number
    if n % 2 == 0:
        return False
    
    # Otherwise go on and  compute r and t
    r, t = 0, n - 1
    while t % 2 == 0:
        r += 1
        t //= 2
        
    # evaluate the test i for k iterations
    for _ in range(k):
        # pick a base a
        a = random.randrange(2, n - 1)
        
        #compute a**t mod n
        x = power(a, t, n)
        if x == 1 or x == n - 1: 
            return True # is there is a non trivial solution 
        else: #else go on
            for _ in range(r - 1):
                x = power(x, 2, n)
                if x == n - 1:
                    break
                else:
                    return False
    return True

######################### SOLOVAY STRASSEN #######################

def solovay_strassen(n, k):
    # corner cases
    if n == 2:
        return True
    if not n & 1:
        return False
    
    # otherwise go on and do the test i for k iterations
    for i in range(k):
        # pick a base a
        a = randrange(2, n - 1)
        # compute Jacobi symbol
        x = legendre(a, n)
        # compute a^((n-1)/2) mod n
        y = power(a, (n - 1) / 2, n)
        
        if (x == 0) or (y != x % n):
            return False  # return Composite
    return True # otherwise return Prime

######################### Evaluation #######################
    
    #evaluate in terms of accuracy
def evaluate_test(method, numbers,composition, k):
    # just to separate the three tests cases
    if "fermat" == method:
        name = "Fermat Probabilistic Test"
        name_test = fermat_probabilistic_test
    elif "m-r"== method:
        name = "Miller Rabin Probabilistic Test"
        name_test = miller_rabin
    elif "s-s" == method:
        name = "solovay Strassen Probabilistic Test"
        name_test = solovay_strassen
    times=[] #list of the running times
    ok2 =[] #list of results
    ok = 0  # to count the correctly classified numbers
    for i in range(len(numbers)): # for every number, compute the test and save how much time it took
        t0= time.time()
        result = name_test(numbers[i], k)
        times.append(float(time.time()-t0))
        if result ==  composition[i]:
            ok = ok +1
            ok2.append(1) # 1 = True (the test is correct) 
        else:
            ok2.append(0) # 0 = False (the test is not correct)
    print('Accuracy of: ', name,   'is ', ok/len(numbers), sep = ' ')
    return times, ok2 # return time needed and results

    #evaluate in terms of accuracy
def eval_factoring(numbers,composition):
    times=[]
    results={number:[] for number in numbers}
    ok2 =[] #list of results
    ok = 0 # to count the correctly classified numbers
    
    # for every number, compute the test and save how much time it took
    for i in range(len(numbers)):
        t0= time.time()
        
        # the multiprocessing is needed to terminate the factoring function in case this takes too much time to terminate and produce a result
        p = multiprocessing.Process(target=factoring, name='factoring', args=(numbers[i],results))
        p.start()
        # Wait 5 seconds for factoring  
        time.sleep(5)
        # Terminate factoring
        p.terminate()
        # Cleanup
        p.join()
        
        times.append(float(time.time()-t0))
        
        if not results[numbers[i]]: #if there isn't any result (i.e. the computation was too expansive and was aborted)
            result = "NAN" #the result is nan
        else:
            result = results[numbers[i]] 
        
        if result ==  composition[i]: #if there is a result that is also correct, save it
            ok = ok +1
            ok2.append(1) # 1 = True (the test is correct)
        elif result == (not composition[i]):
            ok2.append(0) # 0 = False (the test is nit correct)
        else:
            ok2.append(2) # nan (the test has not finished)
            
    print('Accuracy of: Trivial Division is ', ok/len(numbers), sep = ' ')
    print("Computation aborted:",  ok2.count(2), " times")
    return times, ok2 # return time needed and results
    
    
######################### Plotting #######################

            #digits
def barplot(numbers, times, method):
    if "fermat" in method:
        name = "Fermat Probabilistic Test"
    elif "rabin" in method:
        name = "Miller Rabin Probabilistic Test"
    elif "solovay" in method:
        name = "Solovay Strassen Probabilistic Test"
    elif "factoring" in method:
        name = "Factoring Algorithm"
    
    x = [x for x in range(len(numbers))]
    
    # just to adjust the y_axis according to regular numbers and special numbers 
    if len(x)==28: 
        ym = 0.08
    elif len(x)==24:
        ym = 0.04
    else:
        ym = 0.03
    # plot figure   
    plt.figure(figsize = (10,7))
    plt.title(name)
    plt.ylim(ymax = ym, ymin = 0)
    plt.ylabel('Time needed')
    plt.xlabel('# digits')
    plt.xticks(x, numbers, rotation='vertical')
    plt.bar(x, times)

    plt.show()
