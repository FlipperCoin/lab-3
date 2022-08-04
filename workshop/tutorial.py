# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 09:07:21 2022

@author: Lab3
"""

#%%  assignment syntax
a, b = 3, 4

#%% Ex1
# change this code
mystring = "hello"
myfloat = 10.0
myint = 20

# testing code
if mystring == "hello":
    print("String: %s" % mystring)
if isinstance(myfloat, float) and myfloat == 10.0:
    print("Float: %f" % myfloat)
if isinstance(myint, int) and myint == 20:
    print("Integer: %d" % myint)

#%% Ex2 - lists
numbers = []
strings = []
names = ["John", "Eric", "Jessica"]

# write your code here
numbers.append(1)
numbers.append(2)
numbers.append(3)
strings.append("hello")
strings.append('world')
second_name = names[1]


# this code should write out the filled arrays and the second name in the names list (Eric).
print(numbers)
print(strings)
print("The second name on the names list is %s" % second_name)

#%% repeating string
lotsofhellos = "hello" * 10

#%% list operators
print([1,2,3] * 3)
print([1,2] + [3,4,5])

#%% Ex3
x = object()
y = object()

# TODO: change this code
x_list = [x] * 10
y_list = [y] * 10
big_list = x_list + y_list

print("x_list contains %d objects" % len(x_list))
print("y_list contains %d objects" % len(y_list))
print("big_list contains %d objects" % len(big_list))

# testing code
if x_list.count(x) == 10 and y_list.count(y) == 10:
    print("Almost there...")
if big_list.count(x) == 10 and big_list.count(y) == 10:
    print("Great!")
    
 #%% string formatting
name = ["John","Bohn"]
age = 23
print("%s is %d years old." % (name, age))  

#%% Ex4
data = ("John", "Doe", 53.44)
format_string = "Hello %s %s. Your current balance is %s."

print(format_string % data)

#%% in operator in condition
name = "John"
if name in ["John", "Rick"]:
    print("Your name is either John or Rick.")
    
#%% is operator - compare instances
x = [1,2,3]
y = [1,2,3]
print(x == y) # Prints out True
print(x is y) # Prints out False

#%% Ex5
# change this code
number = 100
second_number = 0
first_array = [1,2,2]
second_array = [1,2]

if number > 15:
    print("1")

if first_array:
    print("2")

if len(second_array) == 2:
    print("3")

if len(first_array) + len(second_array) == 5:
    print("4")

if first_array and first_array[0] == 1:
    print("5")

if not second_number:
    print("6")
    
#%% for and ranges
# Prints out the numbers 0,1,2,3,4
for x in range(5):
    print(x)

# Prints out 3,4,5
for x in range(3, 6):
    print(x)

# Prints out 3,5,7
for x in range(3, 8, 2):
    print(x)
    
#%% else for loops - executed on continue or loop condition termination, but not on break
# Prints out 0,1,2,3,4 and then it prints "count value reached 5"

count=0
while(count<5):
    print(count)
    count +=1
else:
    print("count value reached %d" %(count))

# Prints out 1,2,3,4
for i in range(1, 10):
    if(i%5==0):
        break
    print(i)
else:
    print("this is not printed because for loop is terminated because of break but not due to fail in condition")

#%% Ex6
numbers = [
    951, 402, 984, 651, 360, 69, 408, 319, 601, 485, 980, 507, 725, 547, 544,
    615, 83, 165, 141, 501, 263, 617, 865, 575, 219, 390, 984, 592, 236, 105, 942, 941,
    386, 462, 47, 418, 907, 344, 236, 375, 823, 566, 597, 978, 328, 615, 953, 345,
    399, 162, 758, 219, 918, 237, 412, 566, 826, 248, 866, 950, 626, 949, 687, 217,
    815, 67, 104, 58, 512, 24, 892, 894, 767, 553, 81, 379, 843, 831, 445, 742, 717,
    958, 609, 842, 451, 688, 753, 854, 685, 93, 857, 440, 380, 126, 721, 328, 753, 470,
    743, 527
]

# your code goes here
for number in numbers:
    if number % 2 == 0:
        print(number)
    if number == 237:
        break

#%% Ex7
# Modify this function to return a list of strings as defined above
def list_benefits():
    return ["More organized code", "More readable code", "Easier code reuse", "Allowing programmers to share and connect code together"]

# Modify this function to concatenate to each benefit - " is a benefit of functions!"
def build_sentence(benefit):
    return "%s is a benefit of functions!" % benefit

def name_the_benefits_of_functions():
    list_of_benefits = list_benefits()
    for benefit in list_of_benefits:
        print(build_sentence(benefit))

name_the_benefits_of_functions()

#%% creating numpy array from python list
# Create 2 new lists height and weight
height = [1.87,  1.87, 1.82, 1.91, 1.90, 1.85]
weight = [81.65, 97.52, 95.25, 92.98, 86.18, 88.45]

# Import the numpy package as np
import numpy as np

# Create 2 numpy arrays from height and weight
np_height = np.array(height)
np_weight = np.array(weight)

#%% Ex8 - numpy
weight_kg = [81.65, 97.52, 95.25, 92.98, 86.18, 88.45]

import numpy as np

# Create a numpy array np_weight_kg from weight_kg
np_weight_kg = np.array(weight_kg)

# Create np_weight_lbs from np_weight_kg
np_weight_lbs = np_weight_kg * 2.2
# Print out np_weight_lbs
print(np_weight_lbs)

#%% pandas data frame from dict
dict = {"country": ["Brazil", "Russia", "India", "China", "South Africa"],
       "capital": ["Brasilia", "Moscow", "New Dehli", "Beijing", "Pretoria"],
       "area": [8.516, 17.10, 3.286, 9.597, 1.221],
       "population": [200.4, 143.5, 1252, 1357, 52.98] }

import pandas as pd
brics = pd.DataFrame(dict)
print(brics)
print("###")
# set index 
# Set the index for brics
brics.index = ["BR", "RU", "IN", "CH", "SA"]

# Print out brics with new index values
print(brics)

#%% reading csv into pandas frame
# Import pandas as pd
import pandas as pd

# Import the cars.csv data: cars
cars = pd.read_csv('cars.csv')

# Print out cars
print(cars)

#%% indexing pandas dataframe
# Import pandas and cars.csv


print(brics['country'])
print('##')
# Print out country column as Pandas DataFrame
print(brics[['country']])
print('##')
# Print out DataFrame with country and drives_right columns
print(brics[['area', 'country']])

print(brics[1:3])

#%% loc - label indexing, iloc - integer indexing
print(brics.iloc[[1,4]])
print(brics.loc[["BR","SA"]])