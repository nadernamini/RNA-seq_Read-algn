# RNA-Sequencing Read-Alignment Project
## Tiffany Chen and Nader Namini Asl
### CS176-Fall 2018

## Introduction
In this project, we have implemented some basic BWT functions that we covered in class, and then implemented a
simplified version of an aligner for RNA sequencing reads. 


## Rules and Tips

- You can find all required (skeleton) functions/classes/ methods in `project.py`.  You can download all the files as a 
 `.zip` file; _Note_, all the extra files and imported files can be found in the same directory structure as
  `project.py`.
- Python 3 was used to develop the project and only the libraries and modules that are provided by the default 
installation of Python, as well as numpy, have been used.
- The input and outputs of the functions are specific in the docstring for each function.
- A Jupyter Notebook was used for loading in the data and developing/evaluating our alignment algorithm. You can find
a working solution in the Python file `TODO.py`.


## BWT Functions
### Rules and Tips
- The solutions are 0-indexed, and any intervals are [inclusive, exclusive).
- The algorithm was developed with the expectation that it will solely be tested on correctness on small or medium sized
examples (probably length < 10000). Our function do not time out unless it really takes too long to run.
- We assumed an alphabet consisting of `['$', 'A', 'C', 'T', 'G']`, where "$" will always only be the terminator.
- We assumed the `s` parameter in all the function is already terminated by "$".
- For constructing the suffix array, we did not implement the KS algorithm. From past experiences, the Python
implementation of KS was actually slow as hell (probably because of the Python overhead). We used naive sorting on
prefixes instead of radix sort. Our construction of the suffix array is efficient enough to use it for the aligner
part of the project. Also, when we use our function for the aligner, it does not use a radix of over 100 for sorting,
as we assume the memory usage is limited to something reasonable, not to mention that longer radixes take longer
to generate as well.
- Do not delete docstring for exact suffix matches if you want a sanity check of the code.

### Testing
We have provided a simple Python doctest as sanity check for the BWT functions. You can run this by running
```
python -m doctest project.py
```
on the command line.

### Functions
- `get_suffix_array(s)`<br>
Naive implementation of suffix array generation (0-indexed). This code is fast enough so we have enough time in
`Aligner.__init__` (see bottom).

    _Input:_
    - `s`: a string of the alphabet `['A', 'C', 'G', 'T']` already terminated by a unique delimiter '$'
    
    _Output:_
    - list of indices representing the suffix array
    ```
    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    ```
- `get_bwt(s, sa)`<br>
_Input:_
    - `s`: a string terminated by a unique delimiter '$'
    - `sa`: the suffix array of `s`
    
    _Output:_
    - `L`: BWT of `s` as a string

- `get_F(L)`<br>
_Input:_
    - `L = get_bwt(s)`
    
    _Output:_
    - `F`, first column in `Pi_sorted`

- `get_M(F)`<br>
Returns the helper data structure `M` (using the notation from class). `M` is a dictionary that maps character
strings to start indices. i.e. `M[c]` is the first occurrence of `"c"` in `F`.

    If a character `"c"` does not exist in `F`, we set `M[c] = -1`

- `get_occ(L)`<br>
Returns the helper data structure `OCC` (using the notation from class). `OCC` should be a dictionary that maps 
    string character to a list of integers. If `c` is a string character and `i` is an integer, then `OCC[c][i]` gives
    the number of occurrences of character `"c"` in the bwt string up to and including index `i`.
    
- `exact_suffix_matches(p, M, occ)`<br>
Find the positions within the suffix array sa of the longest possible suffix of `p` 
that is a substring of `s` (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    _Input:_
    - `p`: the pattern string
    - `M`, `occ`: buckets and repeats information used by `sp`, `ep`
    
    _Output:_
    - a tuple `(range, length)`
        - `range`: a tuple (start inclusive, end exclusive) of the indices in `sa` that contains
            the longest suffix of `p` as a prefix. `range=None` if no indices matches any suffix of `p`
        - `length`: length of the longest suffix of `p` found in `s`. `length=0` if no indices matches any suffix of `p`

        An example return value would be `((2, 5), 7)`. This means that `p[len(p) - 7 : len(p)]` is
        found in `s` and matches positions `2`, `3`, and `4` in the suffix array.

        ```
        >>> s = 'ACGT' * 10 + '$'
        >>> sa = get_suffix_array(s)
        >>> sa
        [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
        >>> L = get_bwt(s, sa)
        >>> L
        'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
        >>> F = get_F(L)
        >>> F
        '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
        >>> M = get_M(F)
        >>> sorted(M.items())
        [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
        >>> occ = get_occ(L)
        >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
        (True, True, True)
        >>> occ['$']
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        >>> exact_suffix_matches('ACTGA', M, occ)
        ((1, 11), 1)
        >>> exact_suffix_matches('$', M, occ)
        ((0, 1), 1)
        >>> exact_suffix_matches('AA', M, occ)
        ((1, 11), 1)
        ```
        
## RNA Read Alignment