# LinearSankoff

Run LinearSankoff with default setting: 

<code> cat sample.fasta | ./liearsankoff  </code>

Arguments: 

- w: weight on alignment (DEFAULT=0.3)

- b: beam size (DEFAULT=100, 0 for exact search)

- branch: allow to add extra branch in one sequence (DEFAULT=FALSE)


Run LinearSankoff with specific requirements: 

<code> cat sample.fasta | ./liearsankoff -w 0.4 -b 200 --branch </code> 
