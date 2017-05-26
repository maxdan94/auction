# auction
Efficient C implementation of the auction algorithm of Bertsekas for the assignement problem


# Info:

Feel free to use these lines as you wish. This is an efficient C implementation of the Auction algorithm for the assignement problem such as detailed in:  
D.P. Bertsekas, A distributed algorithm for the assignment problem,Laboratory for Information and Decision Systems Working Paper (M.I.T., March 1979).  
It should easily scale to millions of nodes and/or edges if the data is not too adverserial.

# To compile:

"gcc auction.c -O3 -o auction".

# To execute:

"./auction edgelist.txt res.txt eps".
- "edgelist.txt" should contain the edges with integral weights in the bipartite graph: one edge on each line separated by spaces "n1 n2 w", n1 from 0 to n1max and n2 from 0 to n2max.
- "res.txt" contains the results: one edge "n1 n2 w" of the assignement on each line.
- eps is the stepsize, if eps<1/min(n1,n2) then the algorithm is optimal.

Will print some information in the terminal.


# Initial contributors

Maximilien Danisch  
May 2017  
http://bit.ly/maxdan94  
maximilien.danisch@gmail.com
