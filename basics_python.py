#Draft python

#Open command Palette with Ctrl+Shift+P

#Use different environments for different projects (and install packages needed within for each project)
#To create an environment in bash use:
python -m venv .venv
#Activate venv with (also in bash): .venv\Scripts\activate or venv\Scripts\activate or .\venv\Scripts\Activate.ps1
#Python is running when >>> symbol is displayed
#Contr. + S to save
#math test
3+5*4

#Common data types:
#Integer (int) represent positive or negative numbers eg. 3, -321 etc.
#Floating points number (float) represent decimal numbers eg. 3.14, -0.001 etc.
#Character (str) represent text eg. 'Hello', "World" etc.

#Assign value to a variable (variables can't start with a number):
weight_kg = 60
weight_lb = 2.2 * weight_kg

#you can assign multiple values at once:
first, second = 'Grace', 'Hopper'

#To create a string, add single or double quotes.
patient_id = '001'
#To add a prefix:
patient_id = 'inflam_' + patient_id

#display info with:
print(weight_lb)
print(patient_id)

#Calling a function is followed by parenetheses (like in print above)
#(place stuff in quotes to print whatever's written inside)
#for command info use help:
help(print)

#Install and import libraries (installing in venv powershell, importing in python terminal):
pip install numpy
pip install matplotlib
pip install glob
pip install jupyterlab #(open with: jupyter lab)
pip install ipykernel #(to use jupyter notebooks in vscode)
#go to draft_jupyter.ipynb for further info!
import numpy|
import matplotlib.pyplot
import glob

#Lists are inside square brackets (their positions are numbered starting from 0):
odds = [1,3,5,7]
print(odds)
print(odds[2])
#we can also call them in opposite order using minus (here they start from -1)
print(odds[-2])
#To slice a substring of a list, use [start:stop], example:
print(odds[1:3]) #this will print the values at positions 1 and 2 (stop is not included!)

#We can repalce values in lists/strings:
odds[1] = 4
print(odds)
print(len(odds)) #len shows the length of a list

#indices are displayed in [row,column], e.g:
odds[1][3]
#will give the position of a table's cell at row 1 and column 3.
#In a dataframe ":" on its own means "all rows" or "all columns".
#subset command used in dataframes to select specific rows and columns.

#for loop example:
for num in odds:
    print(num)

#if condition example:
num = 37
if num > 100:
    print('greater')
elif num == 0:         #== used for equality, only one "=" is for variables
    print('zero')
else:
    print('not greater')
print('done')

# != (means does not equal)

#
#Define a function using def followed by the name of the function, parameters (in parentheses), and a block of code:
#Defining a function does not run it, similar to a variable.
#Example:
def explicit_fahr_to_celcius(temp):
    converted = ((temp-32)*(5/9))
    return converted

#or more efficiently without creating a variable:
def fahr_to_celcius(temp):
    return (temp-32)*(5/9)

#run example with:
fahr_to_celcius(32)

#Variables within a function only exist within that function and cannot be accessed outside of it.

#Error interpretation:
#Traceback (most recent call last!). Pay attention to the last arrow, the previous ones simply indicate the path of error to previous variables.
#Syntax error (arrow pointing from below) is a simple typo/unknown command.
#NameError occurs when a variable doesn't exist.
#IndexError occurs when a list index does not exist.
#FileError occurs when a file doesn't exist.

#Include assertions in scripts to locate errors when they occur:
assert num > 0, 'number is not positive' #example

#To create a git folder locally (and store history) type this in shell:
git init #This creates the hidden .git folder in the current directory.
#To link it to GitHub repository, we copy the SHH (secure shell protocol) url from there and use it as:
git remote add origin git@github.com:georgerigopoulos98-lab/PCH2a_d90.git
#verify with:
git remote -v
#Then we need to create an SSH key pair:
ssh-keygen -t ed25519 -C "georgerigopoulos98@gmail.com"
#This file is created in C:\Users\georg/ .ssh/id_ed25519
#Passphrase in my notes.
#I can also link R to Git if needed (software carpentry: Version Control with Git 14.)

#Commit only updates my local repository. Push updates the remote repository to update it with the local changes.


#test:
data = numpy.loadtxt(fname='data/inflammation-01.csv', delimiter=',')
image = matplotlib.pyplot.imshow(data)
matplotlib.pyplot.show()


