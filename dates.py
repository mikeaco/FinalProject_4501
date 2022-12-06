import pandas as pd 
import numpy as np 

# dates = []
# with open('new_all.fasta') as topo_file:
#     for line in topo_file:
#         if line[0] == '>':
#             l = line.split()
#             l[0]=l[0].replace('.gb','')
#             dates.append(l)

# #print(dates)

# d = np.array(dates)


# dataframe = pd.DataFrame(d) 
# # dataframe.to_csv(r"dates.csv")

# df = pd.read_csv("dates.csv")
# df.date = df.date.str.replace('JAN', '1') 
# df.date = df.date.str.replace('FEB', '2') 
# df.date = df.date.str.replace('MAR', '3') 
# df.date = df.date.str.replace('APR', '4') 
# df.date = df.date.str.replace('MAY', '5') 
# df.date = df.date.str.replace('JUN', '6') 
# df.date = df.date.str.replace('JUL', '7') 
# df.date = df.date.str.replace('AUG', '8') 
# df.date = df.date.str.replace('SEP', '9') 
# df.date = df.date.str.replace('OCT', '10') 
# df.date = df.date.str.replace('NOV', '11') 
# df.date = df.date.str.replace('DEC', '12')


# df.to_csv(r"dates.csv")

df = pd.read_csv("dates.csv")
s = df.values.tolist()

#print(s)

for i in range(len(s)):
    d = s[i][1].split('-')
    t = d[0]

    d[0] = d[2]
    d[2] = t
    s[i][1] = "-".join(d)
    


df = pd.DataFrame(s) 

df.to_csv(r"dates.csv")
