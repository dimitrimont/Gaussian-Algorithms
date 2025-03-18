
#set matrices for each problem
gjeMTX = [[1, 0, 2, 1],
          [2, -1, 3, -1],
          [4, 1, 8, 2]]

gsMTX = [[1, 0, 2, 1],
         [2, -1, 3, -1],
         [4, 1, 8, 2]]


invertMTX =  [[1, -1, 0],
              [-2, 2, -1],
              [0, 1, -2]]

detMTX =  [[1, -1, 0],
           [-2, 2, -1],
           [0, 1, -2]]

#====================================================================================


def gauss_jordan_elim(mtx1):

    #get row and column length
    r = len(mtx1)
    c = len(mtx1[0])

    #flag
    E = 1

    #go through each column
    for j in range(c - 1):

        #marker in case we need to swap rows
        spot = j

        #get largest mag along diagonal
        largest_col_mag = abs(mtx1[j][j])

        #get largest mag along column and to place on diagonal
        for i in range(j + 1, r):

            #reassign if we find larger mag
            if abs(mtx1[i][j]) > largest_col_mag:
                largest_col_mag = abs(mtx1[i][j])

                #save row to switch to if needed
                spot = i

        #check if elimination is possible
        if largest_col_mag == 0:
            E = 0
            break

        #if the largest is on a row below we swap the rows
        if spot > j:
            mtx1[j], mtx1[spot] = mtx1[spot], mtx1[j]

        #get value along diagonal we will divide each element in the row by
        div = mtx1[j][j]

        #divide each element in the row by the div value
        for x in range(c):
            mtx1[j][x] /= div

        #check if we are not on row we divided by and subtract by the row value times the element postion
        for y in range(r):
            if y != j:
                row_val = mtx1[y][j]
                for z in range(c):
                    #go through each and subtract row value by each element in row
                    mtx1[y][z] -= row_val * mtx1[j][z]

    #return if we have a completed(valid) matrix
    if E != 0:
        return mtx1
    else:
        print("N/A")

#================================================================================================


def gauss_j_inversion(mtx):

    #get row length
    r = len(mtx)

    #create an empty mtx for our identity portion
    ID_mtx = []

    #go through rows
    for i in range(r):
        #this will be for each row in our identity mtx
        identity = []
        #if we are on the diagonal we append 1 else we append 0
        for x in range(r):
            if i == x:
                identity.append(1)
            else:
                identity.append(0)
        #append the identity row to our empty identity matrix
        ID_mtx.append(identity)
    #create a new matrix to store both out identity and our original matrix
    new_mtx = []

    #go through each row and add each line of our original matrix plus our
    #line for identity matrix
    for i in range(r):
        new_mtx.append(mtx[i] + ID_mtx[i])

    #row and column length
    r1 = len(new_mtx)
    c = len(new_mtx[0])

    #flag
    E = 1

    #we dont iterate through the whole matrix
    middle_of_mtx = len(new_mtx[0])//2
    #middle = len(new_mtx[1]) - 1

    #go through each column
    for j in range(middle_of_mtx):

        #marker in case we need to swap rows
        spot = j

        #get largest mag along diagonal
        largest_col_mag = abs(new_mtx[j][j])

        #get largest mag along column and to place on diagonal
        for i in range(j + 1, r1):

            #reassign if we find larger mag
            if abs(new_mtx[i][j]) > largest_col_mag:
                largest_col_mag = abs(new_mtx[i][j])

                #save row to switch to
                spot = i

        #check if elimination is possible
        if largest_col_mag == 0:
            E = 0
            break

        #if the largest is on a row below we swap the rows
        if spot > j:
            new_mtx[j], new_mtx[spot] = new_mtx[spot], new_mtx[j]

        #get pivot we will divide each element in the row by
        div = new_mtx[j][j]

        #divide each element in the row by the div pivot
        for x in range(c):
            new_mtx[j][x] /= div

        #check if we are not on row we divided by and subtract by the row value times the element postion
        for i in range(r1):
            if i != j:
                row_val = new_mtx[i][j]
                for z in range(c):

                    #go through each column and subtract
                    new_mtx[i][z] -= row_val * new_mtx[j][z]

    #return if we have a completed matrix
    if E != 0:
        return new_mtx
    else:
        print("N/A")

#================================================================================================


def gaussian(mtx):

    # row and column length
    r = len(mtx)
    c = len(mtx[0])

    # flag
    E = 1

    #initialize value for counting number row swaps
    swaps = 0

    # go through each column
    for j in range(c - 1):

        # marker in case we need to swap rows
        spot = j

        # get largest mag along diagonal
        largest_col_mag = abs(mtx[j][j])

        #get largest mag along column and to place on diagonal
        for i in range(j + 1, r):

            # reassign if we find larger mag
            if abs(mtx[i][j]) > largest_col_mag:
                largest_col_mag = abs(mtx[i][j])

                # save row to switch to
                spot = i

        # check if elimination is possible
        if largest_col_mag == 0:
            E = 0
            break

        # if the largest is on a row below we swap the rows and add to number of swaps
        if spot > j:
            mtx[j], mtx[spot] = mtx[spot], mtx[j]
            swaps += 1

        for i in range(r):
            if i > j:
                row_val = mtx[i][j]
                for z in range(c):
                    # go through each column and subtract row value divded
                    # by diagonal value times row element
                    mtx[i][z] -= (row_val/mtx[j][j]) * mtx[j][z]

    if E == 1:
        #get rows of mtx(we already have it but just getting it again)
        r1 = len(mtx)

        #create a list to hold the final elements
        answer = []

        #loop rows bottom up
        for j in range(r1 - 1, -1, -1):

            #keep running sum to move to answer list
            total = 0

            #loop through each row
            for i in range(j + 1, r1):

                # add to the total with variable times matching num
                total += mtx[j][i] * answer[i - j - 1]

            #move solved_var to the answer list after dividing it by diagonal
            solved_var = (mtx[j][-1] - total) / mtx[j][j]

            #add to the answers list
            answer.insert(0, solved_var)

        #print answer and return the matrix and number of swaps for determinate function
        print("Gaussian answer: ", answer, "\n")
        return mtx, swaps
    else:
        print("N/A")

#================================================================================================


def mtx_det(mtx):

    #get the matrix and number of swaps by passing mtx to the gaussian fucntion
    new_mtx = gaussian(mtx)[0]
    swaps = gaussian(mtx)[1]

    #print(new_mtx)
    #print(swaps)

    #get columns to go through
    c = len(new_mtx[0])

    #counter for total
    count = 1

    #iterate through each column
    for j in range(c):

        #multiply counter by diagonal values
        count *= new_mtx[j][j]

    #determinate calculation
    det = ((-1)**swaps) * count

    return det

# ================================================================================================

if __name__ == "__main__":
    print("What calculation would you like to perform?\n1. Gauss Jordan Elimination\n"
          "2. Gauss Jordan Inversion\n3. Gaussian Elimination\n4. Matrix Determinate\n5. Exit ")
    user = int(input("Choose a number 1-5: "))

    while user <= 4:
        if user == 1:
            print("Gauss Jorden Elimination Answer: \n")
            for row in gauss_jordan_elim(gjeMTX):
                print(row)

            user = int(input("\nChoose a number 1-5: "))
            #break
        elif user == 2:
            for row in gauss_j_inversion(invertMTX):
                print(row)

            user = int(input("\nChoose a number 1-5: "))
        elif user == 3:
            #print(gaussian(gsMTX))
            gaussian(gsMTX)

            user = int(input("\nChoose a number 1-5: "))
        elif user == 4:
            print("Determinate answer: ", mtx_det(detMTX))

            user = int(input("\nChoose a number 1-5: \n"))
