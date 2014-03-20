#Pallavi Meharia

#Defining the logic variables.
alignment_matrix = []
matrix = []
mylist = []

#Accepting the alignment parameters from the user.
v = input("Please enter the Genome Sequence: ")
v = v.upper()
v_len = len(v)
w = input("Please enter the Genome Sequence: ")
w = w.upper()
w_len = len(w)
match = int(input("Please provide the score against matches: "))
mismatch = int(input("Please provide the mismatch penalty: "))
gap_penalty = int(input("Please provide the gap penalty: "))

#Making the matrix.
def make_list(size):
    mylist = []
    for i in range(size):
        mylist.append(0)
    return mylist
def make_matrix(rows, cols):
    matrix = []
    for i in range(rows):
        matrix.append(make_list(cols))
    return matrix

alignment_matrix = make_matrix(2*w_len+2,v_len+2)

#Provding the column header for the Overlap Alignment Matrix.
x = 0
for i in range(2,v_len+2):
        alignment_matrix[0][i] = v[x]
        x = x + 1
        
#Providing the row header for the Overlap Alignment Matrix.
x = 0
for i in range(2,2*w_len+2,2):
    if i%2 == 0:
        alignment_matrix[i][0] = w[x]
        alignment_matrix[i+1][0] = w[x]
        x = x + 1

##
#Computing the Overlap Alignment Matrix.
x = 0
for j in range(2, 2*w_len+2,2):
    for i in range(2, v_len+2):
        if v[i-2] == w[x]: #match
            alignment_matrix[j][i] = '*'
            alignment_matrix[j+1][i] = alignment_matrix[j-1][i-1] + match
        else: #match not found
            y = alignment_matrix[j-1][i] + gap_penalty
            z = alignment_matrix[j+1][i-1] + gap_penalty
            c = alignment_matrix[j-1][i-1] + mismatch
            mmax = max(y,z,c)
            if y == mmax: #gap 
                alignment_matrix[j][i] = '^'
                alignment_matrix[j+1][i] = y
            if z == mmax: #gap 
                alignment_matrix[j][i] = '<'
                alignment_matrix[j+1][i] = z
            if c == mmax: #mismatch
                alignment_matrix[j][i] = '*'
                alignment_matrix[j+1][i] = c
    x = x + 1
    
#Formatting the Overlap Alignment Matrix.
for i in range(2*w_len+2):
    for j in range(v_len+2):
        if ( j == 0 and i == 1) or ( j == 1 and i == 0):
            alignment_matrix[j][i] = ' '
    if i%2 == 0:
        alignment_matrix[i][0] = ' '

#Displaying the Overlap Alignment Matrix.
print("\nOverlap Alignment Matrix")
print("________________________")
s = [[str(e) for e in row] for row in alignment_matrix]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
print('\n'.join(table))

##
#To find the optimal score value.
#Finding the maximum score along the last column
mylist = []
iRow = len(alignment_matrix) 
jCol = v_len + 1
for i in range(iRow):
    if i%2!=0 and i>1 and alignment_matrix[i-1][jCol] == '*':
        mylist.append(alignment_matrix[i][jCol])
mmax1 = max(mylist)
#Finding the maximum score along the last row
mylist =[]
for j in range(jCol):
    if j>1 and alignment_matrix[iRow-2][j] == '*':
        mylist.append(alignment_matrix[iRow-1][j])
mmax2 = max(mylist)
mmax = max(mmax1,mmax2) #Optimal Score for the given sequences.
print("\nThe Optimal Score is: ",mmax)

##
#Computing the Trace Matrix
A = make_matrix(w_len+1,v_len+1)
print('\n')
x = 0
c = 0
flag_col = 0
flag_row = 0
for i in range(w_len+1):
    for j in range(v_len+1):
        if i == 0 and j!=0:
            A[i][j] = v[x]
            x = x + 1
        if i != 0 and j == 0:
            A[i][j] = w[c]
            c = c + 1
        if i>0 and j>0:
            if alignment_matrix[i*2][j+1] == '*':
                A[i][j] = '*'            
            if alignment_matrix[i*2][j+1] == '*' and (j+1) == (v_len+1) and alignment_matrix[i*2+1][j+1] == mmax:
                A[i][j] = 'M'
                flag_col = 1
            if alignment_matrix[i*2][j+1] == '*' and (i*2+1) == (2*w_len + 1) and alignment_matrix[i*2+1][j+1] == mmax:
                A[i][j] = 'M'
                flag_row = 1
            if alignment_matrix[i*2][j+1] == '^':
                A[i][j] = '|'
            if alignment_matrix[i*2][j+1] == '<':
                A[i][j] = '-'

print("Trace Matrix for Path Estimation. M = Point of Origin")
print("_____________________________________________________\n")
s = [[str(e) for e in row] for row in A]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
print('\n'.join(table))

##
#Path Traceback
print("\nOverlap Alignments with (or without) gaps for Optimal Score computed")
print("____________________________________________________________________\n")
#Finds path if optimal score exists along the last column from the trace matrix.
if flag_col == 1:
    for i in range(w_len,0,-1):
        seq_v = ""
        seq_w = ""
        j = v_len
        if A[i][j] == 'M':
            seq_w = w[i:]
            for x in range(len(seq_w)):
                seq_v = seq_v + '-'
            col = j
            row = i           
            while row != 0:
                if A[row][col] == '*' or A[row][col] == 'M':
                    seq_v = A[0][col] + seq_v
                    seq_w = A[row][0] + seq_w
                    row = row - 1
                    col = col - 1
                if A[row][col] == '-':
                    seq_v = A[0][col] + seq_v
                    seq_w = '-' + seq_w
                    col = col - 1                    
                if A[row][col] == '|':
                    seq_v = '-' + seq_v
                    seq_w = A[row][0] + seq_w
                    row = row - 1
            seq_v = v[:col] + seq_v
            for x in range(col):
                seq_w = '-' + seq_w
            print(seq_v)
            print(seq_w)
            print("\n")

#Finds path if optimal score exists along the last row from the trace matrix.
if flag_col == 1:
    for j in range(v_len,0,-1):
        if j<(v_len):
            seq_v = ""
            seq_w = ""
            i = w_len
            if A[i][j] == 'M':
                seq_v = v[j:]
                for x in range(len(seq_v)):
                    seq_w = seq_w + '-'
                col = j
                row = i
                while col != 0:
                    if A[row][col] == '*' or A[row][col] == 'M':
                        seq_v = A[0][col] + seq_v
                        seq_w = A[row][0] + seq_w
                        row = row - 1
                        col = col - 1
                    if A[row][col] == '-':
                        seq_v = A[0][col] + seq_v
                        seq_w = '-' + seq_w
                        col = col - 1
                    if A[row][col] == '|':
                        seq_v = '-' + seq_v
                        seq_w = A[row][0] + seq_w
                        row = row - 1
                seq_w = w[:row] + seq_w
                for x in range(row):
                    seq_v = '-' + seq_v
                print(seq_v)
                print(seq_w)
                print("\n")
print("End of Alignments")                
