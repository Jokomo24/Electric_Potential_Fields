/*
By means of a matrix or 2D array, this program builds an electric potential charge field between two charged plates 
with sporadically placed source charges with strengths, locations, and dimensions assigned by the user. The program then 
simulates the diffusion (modeled by a differential equation) of electric potential voltages by assigning to each "cell"
the average value of its 3 or 4 neighbors, depending on its location. The program iterates until the field converges to 
a stable state defined by specified difference value between the values of the original matrix and that of the latest 
iteration.

Class: CS 107, Fall 2018
@author Joe Komosa
@version Nov. 9th, 2018
*/
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Point charge struct
typedef struct Charge_struct {
  int      indRow;
  int      indCol;
  double   strength; //fixed potential value associated with the charge
} Charge;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Prints the electric potential matrix
//[in]data[100][100]: Electric potential array
//[in]SizeR: #rows of electric potential array
//[in]SizeC: #cols of electric potential array
//[out/implicit]Electric potential matrix
void printData2D(double data[100][100], int sizeR, int sizeC) {
  int i;
  int j;
  
  for (i = 0; i < sizeR; ++i) {
     for (j = 0; j < sizeC; ++j) {
        printf("%7.1lf ", data[i][j]);
     }
     printf("\n");
  }
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Evolves electric potential array and finds max difference
//[in]data[100][100]: Electric potential array
//[in]SizeR: #rows of electric potential array
//[in]SizeC: #cols of electric potential array
//[in]someCharges[]: User source charge array
//[in]sizeCharge: User source charge array size
//[out]Max difference
//[out/implicit]Evolved electric potential array
double evolveField(double data[100][100],int sizeR, int sizeC, Charge someCharges[], int sizeCharge){
   double dataCopy[100][100];
   int i;
   int j;
   double maxDiff = 0;
//Makes a copy of electric potential array for proper averaging and difference calculation
   for (i = 0; i < sizeR; ++i) {
     for (j = 0; j < sizeC; ++j) {
        dataCopy[i][j] = data[i][j];
     }
   }
//Evolves matrix
   for (i = 1; i < sizeR - 1; ++i) {
      for (j = 1; j < sizeC - 1; ++j) {
      //Inner matrix
         data[i][j] = ((dataCopy[i - 1][j] + dataCopy[i + 1][j] + dataCopy[i][j - 1]+ dataCopy[i][j + 1]) / 4.0);
      //Top row
         data[0][j] = ((dataCopy[1][j] + dataCopy[0][j - 1] + dataCopy[0][j + 1]) / 3.0);
      //Bottom row 
         data[sizeR - 1][j] = ((dataCopy[sizeR - 2][j] + dataCopy[sizeR - 1][j - 1] + dataCopy[sizeR - 1][j + 1]) / 3.0);
      }
   }
//Replaces source charges
   for (i = 0; i < sizeCharge; ++i) {
     data[someCharges[i].indRow][someCharges[i].indCol] = someCharges[i].strength;
   }
//Finds max difference
   for (i = 0; i < sizeR; ++i) {
     for (j = 0; j < sizeC; ++j) {
        if (fabs((dataCopy[i][j] - data[i][j])) > maxDiff) {
           maxDiff = fabs((dataCopy[i][j] - data[i][j]));
        }
      }
   }
   return maxDiff;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Calculates average charge of electric potential array
//[in]data[100][100]: Electric potential array
//[in]SizeR: #rows of electric potential array
//[in]SizeC: #cols of electric potential array
//[out]Average charge value
double calcAvg(double data[100][100],int sizeR, int sizeC) {
  int i;
  int j;
  int count = 0;
  double sumTotal = 0.0;
  double avgValue = 0.0;
  
  for (i = 0; i < sizeR; ++i) {
     for (j = 0; j < sizeC; ++j) {
        sumTotal = sumTotal + data[i][j];
        count++;
     }
  }
   return avgValue = sumTotal / (double) count;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Finds max charge within electric potential matrix
//[in]someCharges[]: Source charge array of type Charge
//[in]sizeCharge: Size of source charge array
//[out]Max source charge of struct type Charge
Charge findMaxCharge(Charge someCharges[], int sizeCharge){
  int i;
  Charge maxCharge;
      maxCharge.indRow = 0;
      maxCharge.indCol = 0;
      maxCharge.strength = 0;
  
  for (i = 0; i < sizeCharge; ++i) {
     if (fabs(someCharges[i].strength) > fabs(maxCharge.strength)) {
        maxCharge.strength = someCharges[i].strength;
        maxCharge.indRow = someCharges[i].indRow;
        maxCharge.indCol = someCharges[i].indCol;
     }
  }
  return maxCharge;
}                                                                                
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//**Additional Created Function**
//Prompts user to provide source charges locations *within matrix* and their strengths
//[in]userCharges: User specified number of source charges
//[in]userChargeArray[]: 1D charge array of type Charge
//[in]sizeR: #rows of electric potential array
//[in]SizeC: #cols of electric potential array
//[out/implicit]Prompt menu and completed user defined 1D charge array
void chargePrompt(int userCharges, Charge userChargeArray[], int sizeR, int sizeC) {
  int i;

  for (i = 0; i < userCharges; ++i) {
     userChargeArray[i].indRow = 0;
     userChargeArray[i].indCol = 0;
     userChargeArray[i].strength = 0.0;
     
     printf("For source #%d\n", i + 1);
     
   //Row index
     printf("Enter row index:\n");
     scanf("%d", &userChargeArray[i].indRow);
     
     while (userChargeArray[i].indRow < 0 || userChargeArray[i].indRow >= sizeR) {
        printf("Index is outside the array, please try again.");
        scanf("%d", &userChargeArray[i].indRow);
     }
   //Column index
     printf("Enter column index:\n");
     scanf("%d", &userChargeArray[i].indCol);
     
     while (userChargeArray[i].indCol < 0 || userChargeArray[i].indCol >= sizeC) {
        printf("Index is outside the array, please try again.");
        scanf("%d", &userChargeArray[i].indCol);
     }
   //Point charge strength
     printf("Enter source strength:\n\n");
     scanf("%lf", &userChargeArray[i].strength);
  }
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Main
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Declared main variables

int main() {
  int      numRows = 0;
  int      numCols = 0;
  double   platePotL = 0.0;
  double   platePotR = 0.0;
  int      numCharges = 0;
  int      i;
  int      j;
  Charge   chargeArray[100];
  Charge   maxCharge;
  double   electPotMatrix[100][100];
  double   maxDiff = 1; //Initialized to 1 so that condition for first iteration loop will be satisfied
  int      iterStep = 0;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Prompts user to provide a matrix size and values for left and right plate potentials for electric potential matrix
  
  printf("Enter # rows:\n");
  scanf("%d", &numRows);
  
  printf("Enter # columns:\n");
  scanf("%d", &numCols);
  
  printf("Enter left plate potential value:\n");
  scanf("%lf", &platePotL);
  
  printf("Enter right plate potential value:\n\n");
  scanf("%lf", &platePotR);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Prompts asking user to provide number of source charges, their locations *within matrix*, and their strengths
  
  printf("Enter # of internal charges (sources of electric potential):\n\n");
  scanf("%d", &numCharges);

//Prompts user until specified number of source charges are assigned their values
  chargePrompt(numCharges, chargeArray, numRows, numCols);
  
//Max charge
  maxCharge = findMaxCharge(chargeArray, numCharges);
  printf("Charge at row = %d and column = %d has potential strength = %.1lf\n\n", maxCharge.indRow, maxCharge.indCol, maxCharge.strength);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Builds initial electric potential matrix

//Initializes all matrix elements to zero
  for (i = 0; i < numRows; ++i) {
     for (j = 0; j < numCols; ++j) {
        electPotMatrix[i][j] = 0.0;
     }
  }
//Sets plate potentials and source charges
  for (i = 0; i < numRows; ++i) {
   //Sets left plate
      electPotMatrix[i][0] = platePotL;
   //Sets right plate
      electPotMatrix[i][numCols - 1] = platePotR;
   //Sets source charges
      if (i < numCharges) {
         electPotMatrix[chargeArray[i].indRow][chargeArray[i].indCol] = chargeArray[i].strength;
      }
  }
//Prints initial electric potential matrix
   printf("Initial Potential Field:\n\n");
   printData2D(electPotMatrix, numRows, numCols);
   printf("\n");
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//Evolves, counts, and prints field at the 1st, 100th, and last iterations
   
   while (maxDiff > .0001) { //continues iterations until max difference no longer satisfies condition
      iterStep++;
      maxDiff = evolveField(electPotMatrix, numRows, numCols, chargeArray, numCharges);
   //1st and 100th iteration
      if (iterStep == 1 || iterStep == 100) {
         printf("Iteration step #%d:\n", iterStep);
         printf("Max Diff = %lf\n", maxDiff);
         printData2D(electPotMatrix, numRows, numCols);
         printf("The average potential is: %lf\n\n", calcAvg(electPotMatrix, numRows, numCols));
      }
   }
//Last iteration
   printf("Iteration step #%d:\n", iterStep);
   printf("Max Diff = %lf\n", maxDiff);
   printData2D(electPotMatrix, numRows, numCols);
   printf("The average potential is: %lf\n", calcAvg(electPotMatrix, numRows, numCols));
   printf("The potential field converged after %d iterations\n\n", iterStep);
   
  return 0;
}