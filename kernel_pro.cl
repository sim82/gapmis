/**
    libgapmis: a library for pairwise sequence aligment with a single gap.
    Copyright (C) 2012 Nikos Alachiotis, Simon Berger, Tomas Flouri, and
    Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

__kernel void gapmis_kernel (__global unsigned int * patsVec, 
			     __global unsigned int * txtsVec, 
			     __global          int * argsVec, 
			     __global          int * txtsLenVec, 
			     __global        float * pensVec,
			     __global          int * hproVec,
			     __global          int * dproVec,
			     __global        float * scrsVec)  
{

	unsigned int groupID = get_group_id(0);
	unsigned int localID = get_local_id(0);

	int i, j, jj, cur_diag_nxt, mis, gap,
	    max_gap = argsVec[2], 
	    blockSize = argsVec[3],
	    maxPatLen = argsVec[4],
	    dproVecGsize = argsVec[5],
	    hproVecGsize = argsVec[6],
	    m = argsVec[groupID + 7],
	    n = txtsLenVec[localID],
	    doffset = 1,gapmismax;


	unsigned int j_min, j_max, abs_ij, patChar, txtChar;

	float temp_score, score = -1000000.0;				
	

int EBLOSUM62_matrix[24][24] =   { { 4, -1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},
  				   {-1, 5,  0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},
  				   {-2, 0, 6,  1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4, -2,-3, 3, 0,-1,-4},
  				   {-2,-2, 1, 6, -3, 0, 2,-1,-1,-3,-4, -1,-3,-3,-1, 0,-1,-4, -3,-3, 4,  1,-1,-4},
  				   {0,-3,-3,-3, 9, -3,-4, -3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
   				   {-1, 1, 0, 0,-3, 5,  2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},
   				   {-1, 0, 0, 2,-4, 2, 5, -2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
  				   {0,-2, 0,-1,-3,-2,-2, 6, -2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},
  				   {-2, 0, 1,-1,-3, 0, 0,-2, 8, -3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},
   				   {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},
  				   {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},
  				   {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5, -1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},
  				   {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,  0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4}, 
  				   {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6, -4,-2,-2, 1, 3,-1,-3,-3,-1,-4},
  				   {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7, -1,-1,-4,-3,-2,-2,-1,-2,-4}, 
   				   {1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},
  				   {0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5, -2,-2, 0,-1,-1, 0,-4},
  				   {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2, 11, 2,-3,-4,-3,-2,-4},
   				   {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7, -1,-3,-2,-1,-4},
  				   {0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},
  				   {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},
  				   {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4}, 
  				   {0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},
  				   {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}};




	if(1 >  m - max_gap)
		score = ( m - 1 ) * pensVec[1] + pensVec[0];

	for( i = 1; i < m; i++ )
	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = max ( 1,  i - max_gap );

		j_max = min ( n,  i + max_gap );

		if(i<= max_gap+1){
			cur_diag_nxt = 0;
		}
		else{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;
		}

		for( j = j_min; j <i; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];			

			gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
			
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;					
								
		}

			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
			

			dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = mis;//gapmismax;//max ( mis, gap );					
								
			

		for( j = i+1; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
			
			gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 
			
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;//max ( mis, gap );				
								
		}		
	}


	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = max ( 1,  i - max_gap );

		j_max = min ( n,  i + max_gap );

		if(i<= max_gap+1){
			cur_diag_nxt = 0;
		}
		else{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;
		}

		for( j = j_min; j < i; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];			
			
			abs_ij = i-j;
				
			gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
	
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;
	
			if(abs_ij<=max_gap)
			{		
				temp_score = gapmismax + ( abs_ij - 1 ) * pensVec[1] + pensVec[0];						
	
				score = max(temp_score,score);	
			}
								
								
		}
	
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
		
			dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;			
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = mis;

		
			if(abs_ij<=max_gap)
			{
				temp_score = mis;		
				
				score = max(temp_score,score);
			}								
	

		for( j = i+1; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];

			abs_ij = j-i;
		
			gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 			
	
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;
		
			if(abs_ij<=max_gap)
			{
				temp_score = gapmismax + ( abs_ij - 1 ) * pensVec[1] + pensVec[0];						

				score = max(temp_score,score);
			}								
		}	
	}
	
	scrsVec[groupID*blockSize+localID]=score;
}


