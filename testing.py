string = "[[-1, -2, 3, 4, 5,6,1,4,5]]"
string = string.replace('[', "")
string = string.replace(']', "")
list = [int(el) for el in string.split(',')]
print(list)
