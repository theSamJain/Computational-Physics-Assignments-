f = eval("lambda x:" + input("Enter the Expression: "))

def secant(f, guess1, guess2, eps): 
    counter = 0       
    x0 = 0 
    run = True

    if (f(guess1) * f(guess2) < 0): 
        
        while (run == True):
            x0 = ((guess1*f(guess2) - guess2*f(guess1))/(f(guess2) - f(guess1)))

            if (f(guess1) * f(x0) < 0):
                guess1 = guess1
                guess2 = x0 
                check0 = f(guess1) * f(x0)

            elif (f(guess2) * f(x0) < 0):
                guess1 = x0
                guess2 = guess2
                check0 = f(guess2) * f(x0)

            counter += 1 

            if (check0 == 0):  
                run = False
                return (x0, counter)
            
            else:
                x_new = ((guess1*f(guess2) - guess2*f(guess1))/(f(guess2) - f(guess1))); 
            
                if (abs((x_new - x0)/x0) <= eps):
                    run = False 
                    return(x0, counter)

                else:
                    while()
                    result = tol()
                    return(result) 
    
    else: 
        print("Can not find a root in the given inteval using Secant Method");  

print("Root of the given equation =", round(x0, 10));  
print("No. of iterations = ", counter); 


guess1 = 0;  
guess2 = 1;  
eps = 0.00001;  
secant(f, guess1, guess2, eps);  