# 18.338: Eigenvalues of random matrices, Fall 2022

## Annoucement

Announcements will be made in piazza most of the times. If you have not received piazza invitation yet, please email sw2030@mit.edu.

   

## Previous versions:
* [2021 Fall 18.338](https://github.com/mitmath/18338/tree/2021-fall)
* [2020 Fall 18.338](https://github.com/mitmath/18338/tree/2020-fall)
* [2019 Fall 18.338](https://github.com/mitmath/18338/tree/2019)
* [2018 Spring 18.338](https://web.mit.edu/18.338/www/2018s/index.html)
* [2016 Spring 18.338](https://web.mit.edu/18.338/www/2016s/index.html)
* [2015 Spring 18.338](https://web.mit.edu/18.338/www/2015s/index.html)
* [2013 Spring 18.338](https://web.mit.edu/18.338/www/2013s/index.html)
* [2012 Spring 18.338](https://web.mit.edu/18.338/www/2012s/index.html)
* [2006 Fall 18.338](https://web.mit.edu/18.338/www/2006f/index.html)
* [2005 Spring 18.325](http://web.mit.edu/18.325/www/)
* [2004 Spring 18.996](https://web.mit.edu/18.338/Spring04/index.html)

## Lecturer: Prof. Alan Edelman

This is the repository for public materials for the MIT course 18.338, *Eigenvalues of random matrices*, for the Fall 2022 semester.

## Location and Time
The class will be held at 2-139 from 3--4:30pm every Monday and Wednesday. 

## Course Description:

We focus on the mathematics of random matrices - from the finite to the infinite, and beyond.

Our emphasis will be on interplay between the varying mathematical tools that have come to play in the modern understanding of random matrix theory. We will also discuss applications of random matrix techniques to problems in engineering and science.

Additional topics will be decided based on the interests of the students. No particular prerequisites are needed though a proficiency in linear algebra and basic probability will be assumed. A familiary with numerical computing languages such as Julia, MATLAB, or Mathematica may be useful .... our primary focus will be Julia and some Mathematica.

This is a graduate course that is intended to be flexible so as to cover the backgrounds of different students. Generally grading will be based on satisfactory completion of problem sets and projects or equivalents.  Homework will be peer graded, and we will look to a rotating student for solutions.

## Homeworks

There will be a few homeworks to start class off including access to the course textbook (being written) where comments are being asked for.
Submit your homework via canvas.mit.edu. If you are not enrolled in 18338 canvas, please email sw2030@mit.edu

|#|Due| Reading| Link |
|-|---|----|------|
| 1 | 9/19 Monday | Ch 5, 9, 10 | [HW 1](https://github.com/mitmath/18338/blob/master/homeworks/fall2022hw1.pdf) |
| 2 | 10/3 Monday | Ch 1, 2, 3 | [HW 2](https://github.com/mitmath/18338/blob/master/homeworks/fall2022hw2.pdf) |
| 3 | 10/17 Monday | Dpp notes | [HW 3](https://github.com/mitmath/18338/blob/master/homeworks/fall2022hw3.pdf) |
| 4 | 10/31 Monday |  | [Project Proposal](https://docs.google.com/document/d/1xW9tE7v1tDmGnN4C1w8cmesiwyeQBaD3HmdLLa6h0ss/edit?usp=sharing) |


## Projects

Please check out [the list of possible class projects](https://github.com/mitmath/18338/blob/master/project_lists.pdf). You are also welcome to work on your own topic. On the second and fourth Fridays, @ 3-4pm, Sungwoo will have office hours for class projects. Students are strongly encouraged to talk to Sungwoo before the end of September, if they are looking for a class project topic. 

There will be a project proposal (maybe a page) due Oct 26 and short 3-minute midterm presentation (project proposal) on 10/31.
The proposal should have references, a short computation or a bit of mathematics or both, something to show that you are on your way.

Please sign up for the final project presentation [here](https://docs.google.com/document/d/11mSnB3_4ad2n4zQ9qlm0Lw5PEyWylJcwKgtAsazh_vM/edit?usp=sharing)


A full project report is due is 12/14.

## Schedule (Tentative)

|#|Day| Date |  Topic | Reading| HW Due |
|-|-|------|------|-----|--|
| 1 | W | 7-Sep |  Hermite, Laguerre and Jacobi ensemble: the ubiquitous triad  | [[Slides]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/HLJslides.pdf)  |  |
| 2 | M | 12-Sep | Semicircle, Quartercircle, Circular and other infinite RMT Laws  |  Ch 3 |    |
| 3 | W | 14-Sep | Random Matrix Decomposition and Finite RMT  | Ch 5  | 
| 4 | M | 19-Sep | Matrix Calculus and Jacobians of Matrix Decompositions  | Ch 9, 10  | HW1  |
| 5 | W | 21-Sep | Determinantal Point Processes (DPP) I  | [[DPPnotes]](https://github.com/mitmath/18338/blob/master/notes/dppnotes.pdf) |  |
| 6 | M | 26-Sep | Determinantal Point Processes (DPP) II  (Sungwoo) |  |  |
| 7 | W | 28-Sep | Algorithms for Sampling DPPs I   |  |  |
| 8 | M | 3-Oct  | Algorithms for Sampling DPPs II  |  | HW2 |
| 9 | W | 5-Oct  | Continuous DPPs (Sungwoo)  | [[Notebook]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/Can%20DPP%20really%20sample%20eigenvalues%3F.ipynb) |  |
|   | M | 10-Oct  | Indigenous People's Day    |  |  |
| 10 | W | 12-Oct | Longest Increasing Sequence (LIS)  |  |  |
| 11 | M | 17-Oct | LIS and RSK algorithm   |  | HW3 |
| 12 | W | 19-Oct | LIS and Group representation   |  |  |
| 13 | M | 24-Oct | LIS and Schur polynomials   |  |  |
| 14 | W | 26-Oct | How many zeros of a random polynomial are real?  | [[Paper]](https://arxiv.org/abs/math/9501224) |  |
| 15 | M | 31-Oct | Project presentation (Midterm)  |  |  |
| 16 | W | 2-Nov |    |  |  |
| 17 | M | 7-Nov | Aztec Diamond and Airy Process    |  |  |
| 18 | W | 9-Nov |     |  |  |
| 19 | M | 14-Nov | Free Probability I   |  |  |
| 20 | W | 16-Nov | Free Probability I   |  |  |
| 21 | M | 21-Nov | Free Probability III : Free Cumulants   |  |  |
|    | W | 23-Nov | *Canceled for Thanksgiving travel*  |  |  |
| 22 | M | 28-Nov | Growth Process I       |  |  |
| 23 | W | 30-Nov | Growth Process II       |  |  |
| 24 | M | 5-Dec  | Growth Process III       |  |  |
| 25 | W | 7-Dec  | Project Presentation I         |  |  |
| 26 | M | 12-Dec  | Project Presentation II       |  |  |
| 27 | W | 14-Dec  | Project Presentation III      |  |  |

## Previous and current Projects 

|Year|Name|Topic|Slides|Report|other materials|
|-|---|-------------|-|-|-|
| 2022 |Xiaomin Li & Yi Tian | Limiting Spectral Distributions of Random Matrices under Finite-Rank Perturbations | | | |
| |Ron Nissim | The KPZ Equation and Fixed Point | | | |
| |Evelyn Ringoot |Largest Singular Values of Bi-Diagonal matrices in Julia | | | |
| | Guanghao Ye & Lichen Zhang |How to Sample Uniform Spanning Trees, DPPs, and NDPPs | | | |
| |Bowen Zhu  |Multivariate Orthogonal Polynomials Symbolically in Julia | | | |
| | | | | | |
| | | | | | |
| | | | | | |
| | | | | | |
| | | | | | |
|2012|Chenhui Hu | Spectral Perturbation of Small-World Networks | | | |
| | Mina Karzand | Random Matrix Theory and Non-coherent MIMO Communications | | | | 
| |Charlotte Kiang | 	RMT Applications to Control Theory | | | |
| |Yi Sun| Virasoro Constraints and the Tracy-Widom Law | | | | 
| |Matt Welborn | Density of states of disordered systems via free addition  | | | | 
| | Helen Xie | RMT Applications to Transport Property | | | | 
| | Yufei Zhao| Spectrual Distributions of Random Graphs | | | |
| |Yi Zeng |North Pole Problem | | | |
 
