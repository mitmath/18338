# 18.338: Eigenvalues of random matrices, Fall 2023

## Announcement

Announcements will be made in piazza most of the times. If you haven't already, please sign up at https://piazza.com/mit/fall2023/18338. 

## Previous versions:
* [2022 Fall 18.338](https://github.com/mitmath/18338/tree/2022-fall)
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

This is the repository for public materials for the MIT course 18.338, *Eigenvalues of random matrices*, for the Fall 2023 semester.

## Location and Time
The class will be held in 2-135 from 3--4:30pm every Monday and Wednesday. 

## Course Description:

We focus on the mathematics of random matrices - from the finite to the infinite, and beyond.

Our emphasis will be on interplay between the varying mathematical tools that have come to play in the modern understanding of random matrix theory. We will also discuss applications of random matrix techniques to problems in engineering and science.

Additional topics will be decided based on the interests of the students. No particular prerequisites are needed though a proficiency in linear algebra and basic probability will be assumed. A familiarity with numerical computing languages such as Julia, MATLAB, or Mathematica may be useful .... our primary focus will be Julia and some Mathematica.

This is a graduate course that is intended to be flexible so as to cover the backgrounds of different students. Generally grading will be based on satisfactory completion of problem sets and projects or equivalents.  Homework will be peer graded, and we will look to a rotating student for solutions.

## Schedule (Tentative)

|#|Day| Date |  Topic | Reading| HW Due |
|-|-|------|------|-----|--|
| 1 | W | 6-Sep |  Hermite, Laguerre and Jacobi ensemble: the ubiquitous triad  | [[Slides]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/HLJslides.pdf)  |  |
| 2 | M | 11-Sep | Semicircle, Quartercircle, Circular and other infinite RMT Laws  |  Ch 3 |    |
| 3 | W | 13-Sep | Random Matrix Decomposition and Finite RMT  | Ch 5  | 
| 4 | M | 18-Sep | Matrix Calculus and Jacobians of Matrix Decompositions  | Ch 9, 10  | HW1  |
| 5 | W | 20-Sep | Determinantal Point Processes (DPP)  | [[DPPnotes]](https://github.com/mitmath/18338/blob/master/notes/dppnotes.pdf) |  |
| 6 | M | 25-Sep | Guest lecture (Alexei Borodin) |  |  |
| 7 | W | 27-Sep | Algorithms for Sampling DPPs   |  |  |
| 8 | M | 2-Oct  | Continuous DPPs | [[Notebook]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/Can%20DPP%20really%20sample%20eigenvalues%3F.ipynb) | HW2 |
| 9 | W | 4-Oct  | Longest Increasing Sequence (LIS)  |  |  |
|   | M | 9-Oct  | Indigenous People's Day    |  |  |
| 10 | W | 11-Oct | LIS and RSK algorithm  |  |  |
| 11 | M | 16-Oct | LIS and Group representation   |  | HW3 |
| 12 | W | 18-Oct | Schur polynomials   |  |  |
| 13 | M | 23-Oct | Jack polynomials I   |  |  |
| 14 | W | 25-Oct | Jack polynomials II  |  |  |
| 15 | M | 30-Oct | Project presentation (Midterm)  |  |  |
| 16 | W | 1-Nov |  Tracy-Widom I  |  |  |
| 17 | M | 6-Nov |  Tracy-Widom II |  |  |
| 18 | W | 8-Nov | Aztec Diamond and Airy Process    |  |  |
| 19 | M | 13-Nov | Growth Process I   |  |  |
| 20 | W | 15-Nov | Growth Process I   |  |  |
| 21 | M | 20-Nov | Growth Process III   |  |  |
|    | W | 22-Nov | *Canceled for Thanksgiving travel*  |  |  |
| 22 | M | 27-Nov | Free Probability I       |  |  |
| 23 | W | 29-Nov | Free Probability II       |  |  |
| 24 | M | 4-Dec  | Free Probability III: Free Cumulants       |  |  |
| 25 | W | 6-Dec  | Project Presentation I         |  |  |
| 26 | M | 11-Dec  | Project Presentation II       |  |  |
| 27 | W | 13-Dec  | Project Presentation III      |  |  |

## Previous Projects 

|Year|Name|Topic|Slides|Report|code,etc..|
|-|---|-------------|-|-|-|
| 2022 |Xiaomin Li & Yi Tian | Limiting Spectral Distributions of Random Matrices under Finite-Rank Perturbations | [link](https://github.com/mitmath/18338/blob/master/projects/2022/XiaominLiandYiTian/xlyt_slides.pdf) |[link](https://github.com/mitmath/18338/blob/master/projects/2022/XiaominLiandYiTian/xlyt_report.pdf) | [link](https://github.com/mitmath/18338/tree/master/projects/2022/XiaominLiandYiTian) |
| |Ron Nissim | The KPZ Equation and Fixed Point | | [link](https://github.com/mitmath/18338/blob/master/projects/2022/rn_report.pdf) | |
| |Evelyn Ringoot |Largest Singular Values of Bi-Diagonal matrices in Julia | [link](https://github.com/mitmath/18338/blob/master/projects/2022/EveylneRingoot/er_slides.pptx) |[link](https://github.com/mitmath/18338/blob/master/projects/2022/EveylneRingoot/er_report.pdf) | [link](https://github.com/mitmath/18338/tree/master/projects/2022/EveylneRingoot) |
| | Guanghao Ye & Lichen Zhang |How to Sample Uniform Spanning Trees, DPPs, and NDPPs | | [link](https://github.com/mitmath/18338/blob/master/projects/2022/gylz_report.pdf) | |
| |Bowen Zhu  |Multivariate Orthogonal Polynomials Symbolically in Julia | |[link](https://github.com/mitmath/18338/blob/master/projects/2022/bz_report.pdf) | |
| | Kyle Fridberg | Random Reflections in R^2 | | [link](https://github.com/mitmath/18338/blob/master/projects/2022/kf_report.pdf) | |
| |Luke Robitaille | Moments of Wachter Law | | [link](https://github.com/mitmath/18338/blob/master/projects/2022/lr_report.pdf) | |
| 2021  | Aleksandr Zimin | The Weingartens function for beta=1,2,4 and its possible generalizations using Jack polynomials| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Aleksandr%20Zimin/Jack_polynomials.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Aleksandr%20Zimin/Jack%20polynomials%20project.ipynb)| |
| | Aviva Englander | Maximum Eigenvalues in Brownian Motion and Their Correlation with the Airy Process | | [link](https://github.com/mitmath/18338/blob/master/projects/2021/Aviva%20Englander/Final_Report18.338.pdf)| [link](https://github.com/mitmath/18338/tree/master/projects/2021/Aviva%20Englander)|
| | Daniel Pickard |Sampling Continuous Determinantal Point Processes with ApproxFun |[link](https://github.com/mitmath/18338/blob/master/projects/2021/Daniel%20Pickard/presentation.pdf) |[link](https://github.com/mitmath/18338/blob/master/projects/2021/Daniel%20Pickard/daniel_pickard_term_project.pdf) |[link](https://github.com/mitmath/18338/tree/master/projects/2021/Daniel%20Pickard) |
| |Hussein Fellahi | Random Matrix Filtering| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Hussein%20Fellahi/18338_Presentation-2.pdf) |[link](https://github.com/mitmath/18338/blob/master/projects/2021/Hussein%20Fellahi/18_338_Final_Project.pdf) | [link](https://github.com/mitmath/18338/tree/master/projects/2021/Hussein%20Fellahi)|
| |Jiahai Feng |Hypothesis testing in high dimensions |[link](https://github.com/mitmath/18338/blob/master/projects/2021/Jiahai%20Feng/slidespdf.pdf) |[link](https://github.com/mitmath/18338/blob/master/projects/2021/Jiahai%20Feng/report.pdf) | |
| | Josefina Menendez| What is the expected number of points drawn from a Determinantal Point Process defined by a Wishart kernel? | | |[link](https://github.com/mitmath/18338/tree/master/projects/2021/Josefina%20Menendez) |
| | Madhav Sankaranarayanan| Determinantal Point Processes and Growth Models| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Madhav%20Sankaranarayanan/18.338%20Project%20Presentation.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Madhav%20Sankaranarayanan/18.338%20Project%20Report-1.pdf)|[link](https://github.com/mitmath/18338/tree/master/projects/2021/Madhav%20Sankaranarayanan) |
| | Saaketh Vedantam| Analyzing Higher Order Effects on Eigenvalues| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Saaketh%20Vedantam/Presentation.pdf)|[link](https://github.com/mitmath/18338/blob/master/projects/2021/Saaketh%20Vedantam/Final%20Report.pdf) | |
| |Theo Diamandis | Randomization to Speed up Convex Optimization| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Theo%20Diamandis/presentation.jl)| [link](https://github.com/mitmath/18338/blob/master/projects/2021/Theo%20Diamandis/tdiamand-report.pdf)| |
| 2020 |Max Li |Exploring Densities of Gaussian Quadratic Forms | [link](https://github.com/mitmath/18338/blob/master/projects/2020/ML_slides.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2020/ML_report.pdf)| |
| |Poorya Habibzadeh |Deriving a closed form for the Cauchy transform of two laws | |[link](https://github.com/mitmath/18338/blob/master/projects/2020/PH_report.pdf) | |
| |Shawn Im |Determinantal Point Processes and β-ensembles | |[link](https://github.com/mitmath/18338/blob/master/projects/2020/SI_report.pdf) | |
| |Chun-Hei Lam |Computation of Equilibrium Measure | [link](https://github.com/mitmath/18338/blob/master/projects/2020/SL_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2020/SL_report.pdf) | |
| | Tony Tohme| Simplexmethod and random matrices| [link](https://github.com/mitmath/18338/blob/master/projects/2020/TT_slides.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2020/TT_report.pdf)| |
| 2016 |Antoni Musolas |Differential geometrical approach to covariance estimation | | | |
| |Manishika Agaskar |	Diffraction gain of free-space optical communications in atmospheric turbulence | | | |
| |Hong Hu | Spectral Initialization and its Performance Analysis| | | |
| |Anuran Makur	|Maximal Correlation Functions: Hermite, Laguerre, and Jacobi  | | | |
| | John Urschel |On the Minimal Eigenpair of Erdos-Renyi Graphs | | | |
| |John Napp |RMT and the complexity of linear optics | | | |
| |Alex Wein |Random Matrix Contiguity | | | |
| |Ravi Bajaj |Central Limit Theorem for log-Determinant of Wigner Matrices | | | |
| |Brandon Tran |Application of phase transitions in spiked covariance matrices to synchronization problems | | | |
|2012|Chenhui Hu | Spectral Perturbation of Small-World Networks | | | |
| | Mina Karzand | Random Matrix Theory and Non-coherent MIMO Communications | | | | 
| |Charlotte Kiang | 	RMT Applications to Control Theory | | | |
| |Yi Sun| Virasoro Constraints and the Tracy-Widom Law | | | | 
| |Matt Welborn | Density of states of disordered systems via free addition  | | | | 
| | Helen Xie | RMT Applications to Transport Property | | | | 
| | Yufei Zhao| Spectral Distributions of Random Graphs | | | |
| |Yi Zeng |North Pole Problem | | | |
 
