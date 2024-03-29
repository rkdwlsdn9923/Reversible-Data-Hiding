#include <stdio.h>
#include <math.h>

int main(void)
{
	int coefficient[37][4] = // QDCT coefficient - alter for freq of occurrence
	{ {0, 0, 0, 0},
	{0, -1, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, -1, 0},
	{0, 0, 0, 1},
	{0, 0, 0, -1},
	{1, 0, 0, 0},
	{-1, 0, 0, 0},
	{0, -1, 1, 0},
	{0, 1, -1, 0},
	{0, -1, -1, 0},
	{0, 1, 1, 0},
	{0, 0, -1, 1},
	{0, 0, 1, -1},
	{0, 0, -1, -1},
	{0, 0, 1, 1},
	{-1, 1, 0, 0},
	{1, -1, 0, 0},
	{1, 1, 0, 0},
	{-1, -1, 0, 0},
	{0, 1, 0, 1},
	{0, -1, 0, -1},
	{0, -1, 0, 1},
	{0, 1, 0, -1},
	{1, 0, 1, 0},
	{-1, 0, -1, 0},
	{-1, 0, 1, 0},
	{1, 0, -1, 0},
	{1, 0, 0, -1},
	{1, 0, 0, 1},
	{-1, 0, 0, 1},
	{-1, 0, 0, -1},
	{0, -1, 1, -1},
	{0, -1, -1, -1},
	{0, 1, -1, 1},
	{0, 1, 1, 1}, };

	int secret[16][4] = // secret data bitstream - extended
	{ {0, 0, 0, 0},
	{0, 0, 0, 1},
	{0, 0, 1, 0},
	{0, 0, 1, 1},
	{0, 1, 0, 0},
	{0, 1, 0, 1},
	{0, 1, 1, 0},
	{0, 1, 1, 1},
	{1, 0, 0, 0},
	{1, 0, 0, 1},
	//{1, 0, 1},
	{1, 0, 1, 0},
	{1, 0, 1, 1},
	{1, 1, 0, 0},
	{1, 1, 0, 1},
	//{1, 1, 1}
	{1, 1, 1, 0},
	{1, 1, 1, 1}, };

	int save[37][16][4] = { 0, }; // save output initialize
	for (int i = 0; i < 37; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				save[i][j][k] = coefficient[i][k]; // copy QDCT
			}
		}
	}

	//compensation embedding
	for (int i = 0; i < 37; i++) // i's coefficient test
	{
		for (int j = 0; j < 16; j++) // secret message
		{
			for (int k = 0; k < 4; k++) // r7, r8, r9, r10
			{
				int dh = (int)ceil(((double)save[i][j][(k % 4)] - (double)save[i][j][(k + 1) % 4]) / (double)2);
				if (dh == 0)
				{
					save[i][j][k % 4] += secret[j][k];
					save[i][j][(k + 1) % 4] -= secret[j][k];
				}
				else if (dh > 0)
				{
					save[i][j][k % 4] += 1;
					save[i][j][(k + 1) % 4] -= 1;
				}
			}
		}
	}

	# define score_test 1 // score search
	int score_1[81][5] = // coefficient score from excel // FPP(0.9) + PSNR(0.1)
	{
	{0,0,0,0,1},
	{1,0,0,0,4},
	{-1,0,0,0,5},
	{0,1,0,0,6},
	{0,-1,0,0,7},
	{0,0,1,0,8},
	{0,0,-1,0,8},
	{1,-1,0,0,9},
	{-1,1,0,0,9},
	{0,0,0,1,10},
	{0,0,0,-1,10},
	{0,1,-1,0,13},
	{0,-1,1,0,13},
	{1,1,0,0,14},
	{-1,-1,0,0,15},
	{1,0,-1,0,16},
	{-1,0,1,0,17},
	{0,1,1,0,18},
	{0,-1,-1,0,18},
	{1,0,0,-1,20},
	{-1,0,0,1,20},
	{1,0,1,0,21},
	{-1,0,-1,0,22},
	{0,-1,0,1,24},
	{0,1,0,-1,25},
	{1,0,0,1,25},
	{-1,0,0,-1,26},
	{0,0,-1,1,27},
	{0,0,1,-1,27},
	{0,1,0,1,29},
	{0,-1,0,-1,30},
	{0,0,1,1,32},
	{0,0,-1,-1,32},
	{1,-1,1,0,33},
	{1,1,-1,0,33},
	{-1,1,-1,0,34},
	{-1,1,1,0,34},
	{-1,-1,1,0,34},
	{1,-1,-1,0,34},
	{1,1,1,0,37},
	{-1,-1,-1,0,37},
	{1,-1,0,1,40},
	{-1,1,0,-1,41},
	{1,1,0,-1,41},
	{-1,1,0,1,41},
	{-1,-1,0,1,42},
	{1,-1,0,-1,42},
	{1,1,0,1,44},
	{-1,-1,0,-1,45},
	{0,-1,1,1,48},
	{0,1,-1,1,49},
	{0,1,-1,-1,49},
	{0,1,1,-1,49},
	{0,-1,1,-1,49},
	{0,-1,-1,1,49},
	{0,1,1,1,52},
	{0,-1,-1,-1,52},
	{1,0,-1,1,54},
	{-1,0,1,-1,55},
	{1,0,1,-1,55},
	{-1,0,1,1,55},
	{1,0,-1,-1,55},
	{-1,0,-1,1,55},
	{1,0,1,1,58},
	{-1,0,-1,-1,59},
	{1,-1,1,-1,60},
	{1,1,-1,-1,60},
	{1,-1,-1,1,61},
	{1,1,1,-1,65},
	{1,-1,1,1,65},
	{1,1,-1,1,65},
	{1,-1,-1,-1,65},
	{1,1,1,1,66},
	{-1,1,-1,1,67},
	{-1,-1,1,1,67},
	{-1,1,1,-1,69},
	{-1,1,1,1,72},
	{-1,-1,-1,1,72},
	{-1,1,-1,-1,73},
	{-1,-1,1,-1,73},
	{-1,-1,-1,-1,74},
	};
	int score_2[81][5] = // coefficient score from excel // FPP(0.95) + PSNR(0.05)
	{
	{0,0,0,0,1},
	{1,0,0,0,3},
	{-1,0,0,0,3},
	{0,1,0,0,5},
	{0,-1,0,0,5},
	{0,0,1,0,6},
	{0,0,-1,0,7},
	{0,0,0,1,9},
	{0,0,0,-1,9},
	{1,-1,0,0,9},
	{-1,1,0,0,9},
	{1,1,0,0,12},
	{-1,-1,0,0,12},
	{0,1,-1,0,13},
	{0,-1,1,0,13},
	{0,1,1,0,16},
	{0,-1,-1,0,16},
	{1,0,-1,0,17},
	{-1,0,1,0,17},
	{1,0,1,0,19},
	{-1,0,-1,0,19},
	{1,0,0,-1,21},
	{-1,0,0,1,21},
	{1,0,0,1,23},
	{-1,0,0,-1,24},
	{0,-1,0,1,25},
	{0,1,0,-1,25},
	{0,1,0,1,28},
	{0,-1,0,-1,28},
	{0,0,-1,1,28},
	{0,0,1,-1,28},
	{0,0,1,1,31},
	{0,0,-1,-1,31},
	{1,-1,1,0,33},
	{1,1,-1,0,33},
	{-1,1,-1,0,34},
	{-1,1,1,0,34},
	{-1,-1,1,0,34},
	{1,-1,-1,0,34},
	{1,1,1,0,35},
	{-1,-1,-1,0,35},
	{1,-1,0,1,41},
	{-1,1,0,-1,41},
	{1,1,0,-1,41},
	{-1,1,0,1,41},
	{-1,-1,0,1,42},
	{1,-1,0,-1,42},
	{1,1,0,1,43},
	{-1,-1,0,-1,43},
	{0,-1,1,1,49},
	{0,1,-1,1,49},
	{0,1,-1,-1,49},
	{0,1,1,-1,50},
	{0,-1,1,-1,50},
	{0,-1,-1,1,50},
	{0,1,1,1,51},
	{0,-1,-1,-1,51},
	{1,0,-1,1,56},
	{-1,0,1,-1,56},
	{1,0,1,-1,56},
	{-1,0,1,1,56},
	{1,0,-1,-1,56},
	{-1,0,-1,1,56},
	{1,0,1,1,58},
	{-1,0,-1,-1,58},
	{1,-1,1,-1,63},
	{1,1,-1,-1,63},
	{1,-1,-1,1,63},
	{1,1,1,-1,65},
	{1,-1,1,1,65},
	{1,1,-1,1,65},
	{1,-1,-1,-1,65},
	{1,1,1,1,66},
	{-1,1,-1,1,70},
	{-1,-1,1,1,71},
	{-1,1,1,-1,71},
	{-1,1,1,1,73},
	{-1,-1,-1,1,73},
	{-1,1,-1,-1,73},
	{-1,-1,1,-1,73},
	{-1,-1,-1,-1,74},
	};

	// output display
	for (int i = 0; i < 37; i++) // i's coefficient test
	{
		for (int j = 0; j < 16; j++) // secret message
		{
			for (int k = 0; k < 4; k++) // r7, r8, r9, r10
			{
				if (k == 0 && j == 0) printf("\t  [%d %d %d %d] \n", coefficient[i][0], coefficient[i][1], coefficient[i][2], coefficient[i][3]); // title display				
				printf("%d ", save[i][j][k]); // main display
			}
#if score_test // score search
			int temp_score;
			for (int s = 0; s < 81; s++)
			{
				if (score[s][0] == save[i][j][0] && score[s][1] == save[i][j][1] && score[s][2] == save[i][j][2] && score[s][3] == save[i][j][3])
				{
					temp_score = score[s][4];
					break;
				}
				else temp_score = 0;
			}
			if (temp_score) printf("\t%d", temp_score); // score display
#endif
			printf("\n");
		}
		printf("\n");
	}

	return 0;
}