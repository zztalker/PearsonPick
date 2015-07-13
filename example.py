def fun(x):
    return x*x

x = [1,2,3,4,5]
y = []
for xi in x:
    y += [fun(xi)]

print(x,y)
