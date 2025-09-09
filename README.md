# 18.338: Eigenvalues of random matrices, Fall 2025


## Announcement
We will only use canvas for homework submission.  This github repo is the official course home page.
Announcements will be made in piazza most of the time. If you haven't already, please sign up at https://piazza.com/mit/fall2025/18338. 

## Previous versions:
* [2024 Fall 18.338](https://github.com/mitmath/18338/tree/2024-fall)
* [2023 Fall 18.338](https://github.com/mitmath/18338/tree/2023-fall)
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

## TA: Nicholas West

This is the repository for public materials for the MIT course 18.338, *Eigenvalues of random matrices*, for the Fall 2025 semester.

## Location and Time
The class will be held in 4-153 (in 2025) from 3--4:30pm every Monday and Wednesday. 

## Course Description:

We focus on the mathematics of random matrices - from the finite to the infinite, and beyond.
I'd like this year to cover the basic RMT just quickly and slowly enough to get to 
exciting developments in DPPs, growth processes, aztec diamonds, and some combinatorial structures.

Our emphasis will be on interplay between the varying mathematical tools that have come to play in the modern understanding of random matrix theory. We will also discuss applications of random matrix techniques to problems in engineering and science. In 2024 we may focus on orthogonal polynomial theory. 

Additional topics will be decided based on the interests of the students. No particular prerequisites are needed though a proficiency in linear algebra and basic probability will be assumed. A familiarity with numerical computing languages such as Julia, MATLAB, or Mathematica may be useful .... our primary focus will be Julia and some Mathematica.

This is a graduate course that is intended to be flexible so as to cover the backgrounds of different students. Generally grading will be based on satisfactory completion of problem sets and projects or equivalents.  Homework may be peer graded, and we may look to a rotating student for solutions.  Homework may be every week or two.

## Homeworks 

There will be a few homeworks to start class off including access to the course textbook (being written) where comments are being asked for.
Submit your homework via canvas.mit.edu. If you are not enrolled in 18338 canvas, please email the TA.

|#|Due| Reading| Link |
|-|---|----|------|
| 1| 9/17 | 5,10,11| [HW1](https://github.com/mitmath/18338/blob/master/homeworks/fall2025hw1.pdf)|


## Research Questions
|Topics|
|-|
|[Combinatorics of Wachter's Law](https://docs.google.com/presentation/d/1CelJYt1GOMG8_HlvWqKK7__zFOSXWUGxsvOPwAxqTGo/edit?usp=sharing)|
| [pdf listing from 2023](https://github.com/mitmath/18338/blob/master/project_lists/project_lists.pdf) |
MOPS in Julia
Lanczos on Multivariate Orthogonal Polynomials
QR Multivariate implementation (in Julia, Symbolically if possible)
Check Hermite, Laguerre, Jacobi Eigenvalue (of Laplace-Beltrami) - check it in Julia
second derivative of det(moment) and painleve
Simple derivation paper applied to Laguerre and Jacobi
See if there’s any easier derivation of TW starting from their R, S recursion
Simple TW derivation for GOE/GSE
Simple TW derivation for LUE/JUE

# Monday & Wednesday Schedule (Tentative)
**September 3 - December 10, 2025**

| # | Day | Date | Topic | Reading | HW Due |
|---|-----|------|-------|---------|--------|
| 1 | M | 09/03 | Hermite, Laguerre and Jacobi ensemble: the ubiquitous triad  | [[Slides]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/mit_18338_lecture1.pdf)  [[Growth Notebook]](https://mitmath.github.io/18338/NotebooksReferencedBook/RandomGrowth.html) | |
| 2 | W | 09/05 | Semicircle &  TracyWidom  Laws (Fredholm Det, Diff Eq, & Wigner's Proof) | 1,3 | |
| 3 | M | 09/08 | Random Matrix Decomposition and Finite RMT | | |
| 4 | W | 09/10 | Orthogonal Polynomial Intro  | | |
| 5 | M | 09/15 | Introduction to Discrete DPP   | | |
| 6 | W | 09/17 | DPPs in Random Matrix Theory  | | |
| 7 | M | 09/22 |  Longest Increasing Subsequence | | |
| 8 | W | 09/24 | Guest Lecturer Alexei Borodin | | |
| 9 | M | 09/29 | Growth Processes, Aztec Diamonds and Tracy Widom     | | |
| 10 | W | 10/01 | Brownian Motion and Airy Processes | | |
| 11 | M | 10/06 | catch up | | |
| 12 | W | 10/08 | | | |
| - | M | 10/13 | Indigenous Peoples Day | | |
| 13 | W | 10/15 | | | |
| 14 | M | 10/20 | | | |
| 15 | W | 10/22 | | | |
| 16 | M | 10/27 | | | |
| 17 | W | 10/29 | | | |
| 18 | M | 11/03 | | | |
| 19 | W | 11/05 | | | |
| - | M | 11/10 | Student Holiday | | |
| 20 | W | 11/12 | | | |
| 21 | M | 11/17 | | | |
| 22 | W | 11/19 | | | |
| 23 | M | 11/24 | | | |
| 24 | W | 11/26 | | | |
| 25 | M | 12/01 | | | |
| 26 | W | 12/03 | | | |
| 27 | M | 12/08 | | | |
| 28 | W | 12/10 | | | |

## Some old Schedules 

|#|Day| Date |  Topic | Reading| HW Due |
|-|-|------|------|-----|--|
| 1 | W | 4-Sep |  Hermite, Laguerre and Jacobi ensemble: the ubiquitous triad  | [[Slides]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/mit_18338_lecture1.pdf)  |  |
| 2 | M | 9-Sep | Semicircle, Quartercircle, Circular and other infinite RMT Laws  |  Ch 3  [Research: Wachter's Law Combinatorics](https://docs.google.com/presentation/d/1CelJYt1GOMG8_HlvWqKK7__zFOSXWUGxsvOPwAxqTGo/edit?usp=sharing) |    |
| 3 | W | 11-Sep | Random Matrix Decomposition and Finite RMT  | Ch 5  | 
| 4 | M | 16-Sep | Matrix Calculus and Jacobians of Matrix Decompositions | Ch 10, 11 | [HW 1](https://github.com/mitmath/18338/blob/master/homeworks/fall2024hw1.pdf) Due|
| 5 | W | 18-Sep | Matrix Calculus and Differential Forms   |  Ch 11 |  |
| 6 | M | 23-Sep | Orthogonal Polynomial Intro  | Project Idea: MOPS in Julia |  |
| 7 | W | 25-Sep |Multivariate Orthogonal Polynomial Theory   |  |  |
| 8 | M | 30-Sep |  Univariate Orthogonal Polynomial Theory  |  |  |
| 9 | W | 2-Oct  |   Orthogonal Polynomials and RMT  | | |
|10|  W | 7-Oct | Introduction to Discrete DPP |||
| 11 | W | 9-Oct | DPPs in Random Matrix Theory |||
|    | M | 14-Oct  | Indigenous People's Day    |  |  |
| 12 | W | 16-Oct | Possible Class Projects  | | |
| 13 | M | 21-Oct | Computational Experiments with DPP | | |
| 14 | W | 23-Oct | Algorithms for DPPs | | |
| 15 | M | 28-Oct | Projection DPP Algorithm / Kesten-McKay Law | | |
| 16 | W | 30-Oct | Growth Processes, Aztec Diamonds and Tracy Widom | | |
| 17 | M | 4-Nov |  Longest Increasing Subsequence |  |  |
| 18 | W | 6-Nov | Brownian Motion | |  |
|    | M | 11-Nov | Veterans Day  |  |  |
| 19 | W | 13-Nov |  Project Plan Presentations |  |  |
| 20 | M | 18-Nov | Free Probability Intro |  |  |
| 21 | W | 20-Nov | Free Probability R-Transform |  |  |
| 22 | M | 25-Nov | How many roots are real? |  |  |
|    | W | 27-Nov | *Canceled for Thanksgiving travel*  |  |  |
| 23 | M | 2-Dec |
| 24 | W | 4-Dec  | Project Presentation I         |  |  |
| 25 | M | 9-Dec  | Project Presentation II       |  |  |
| 26 | W | 11-Dec  | Project Presentation III      |  |  |
| | |  We are here ↑ ↑ Below is tentative  and maybe should be ignored| ||  |
| 6 | M | 23-Sep | Multivariate Orthogonal Polynomials   |   [Aztec Notes](https://github.com/mitmath/18338/blob/master/notes/borodin_aztec.pdf)  | [HW 2](https://github.com/mitmath/18338/blob/master/homeworks/fall2024hw2.pdf) Due |
| 7 | W | 25-Sep | Algorithms for Sampling DPPs   |  |  |
| 8 | M | 30-Sep  |DPPs | [[Notebook]](https://github.com/mitmath/18338/blob/master/JuliaNotebooks/Can%20DPP%20really%20sample%20eigenvalues%3F.ipynb)  [[DPPnotes]](https://github.com/mitmath/18338/blob/master/notes/dppnotes.pdf) |  |
| 9 | W | 2-Oct  | Longest Increasing Sequence (LIS)  |  |  |
| 10 | M | 7-Oct | LIS and RSK algorithm  |  |  |
| 11 | W | 9-Oct | LIS and Group representation   |  |  |
|    | M | 14-Oct  | Indigenous People's Day    |  |  |
| 12 | W | 16-Oct | Schur polynomials   |  | [HW 3](https://github.com/mitmath/18338/blob/master/homeworks/fall2023hw3.pdf) Due  |
| 13 | M | 21-Oct | Jack polynomials I   |  |  |
| 14 | W | 23-Oct | Jack polynomials II  |  |  |
| 15 | M | 28-Oct | Project presentation (Midterm)  |  |  |
| 16 | W | 30-Oct |  Tracy-Widom I  |  |  |
| 17 | M | 4-Nov |  Tracy-Widom II |  |  |
| 18 | W | 6-Nov | Aztec Diamond and Airy Process    |  |  |
|    | M | 11-Nov | Veterans Day  |  |  |
| 19 | W | 13-Nov | Growth Process I   |  |  |
| 20 | M | 18-Nov | Growth Process III   |  |  |
| 21 | W | 20-Nov | Free Probability I       |  |  |
| 22 | M | 25-Nov | Free Probability II       |  |  |
|    | W | 27-Nov | *Canceled for Thanksgiving travel*  |  |  |
| 23   | M | 2-Dec  | Free Probability III: Free Cumulants       |  |  |
| 24 | W | 4-Dec  | Project Presentation I         |  |  |
| 25 | M | 9-Dec  | Project Presentation II       |  |  |
| 26 | W | 11-Dec  | Project Presentation III      |  |  |

## Previous Projects 

|Year|Name|Topic|Slides|Report|code,etc..|
|-|---|-------------|-|-|-|
| 2024 | Cecelia Chen | Computing Multivariate Orthogonal Bases ||[link](https://github.com/mitmath/18338/blob/master/projects/2024/chencecilia_100854_4907163_18_338_project%20(1).pdf)|[link](https://github.com/mitmath/18338/blob/master/projects/2024/chencecilia_100854_4907164_MultivariateLanczos.jl%20copy.zip)
| | Vaibhav Dixit | Evolution of Eigenvalue Spectrum of Fully Connected (Dense) Layers in DNNs || [link](https://github.com/mitmath/18338/blob/master/projects/2024/dixitvaibhav_153380_4902546_blan%20(1).pdf)|
| | Jennifer Hritz | Constructing approximately Haar-random unitaries from GUEs| [link](https://github.com/mitmath/18338/blob/master/projects/2024/hritzjennifer_161909_4911299_final_project_presentation_hritz.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/hritzjennifer_161909_4911297_final_project_hritz.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/hritzjennifer_161909_4911298_project_plots.jl)
| | Joonsoo Lee | Derivatives of Randomized Matrix Approximation Algorithms | | [link](https://github.com/mitmath/18338/blob/master/projects/2024/leejoonsoo_174493_4904826_18_338_Project-1.pdf)|
| | Donald Stralkus | Universality in Inner-Product Random Geometric Graphs | | [link](https://github.com/mitmath/18338/blob/master/projects/2024/stralkusdonald_119864_4906139_18338_Final_Project.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/stralkusdonald_119864_4906140_graphSpectra.jl) |
| | Nicholas West | Fast eigendecomposition of unitary upper Hessenberg matrices | [link](https://github.com/mitmath/18338/blob/master/projects/2024/RMT_Project_Presentation.pdf) |[link](https://github.com/mitmath/18338/blob/master/projects/2024/RMT_Report_WEST.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/FastUnitaryEigenvalues.zip)
| | Joanna Zou | k-DPPs: Fixed-Size Determinantal Point Processes for Diversity-Based Subsampling | [link](https://github.com/mitmath/18338/blob/master/projects/2024/zoujoanna_129477_4904902_18338_project_presentation.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/zoujoanna_129477_4904901_18_338_Final_Project_Report.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/zoujoanna_129477_4904903_final_project.ipynb)
| | Alec Zhu |  Spacing of Riemann Zeta Zeros in Julia | [link](https://github.com/mitmath/18338/blob/master/projects/2024/Riemann%20Zeta%20Zeros%20in%20Julia.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2024/zhualec_101690_4905049_submission.zip)|
| 2023 |Gaurav Arya| Differentiable Determinental Point Processes |  |[link](https://github.com/mitmath/18338/blob/master/projects/2023/aryagaurav_LATE_101101_3887118_main.pdf)| 
| |Andrey Bryutkin & Diego Chavez| Painlevé Systems and Eigenvalue Distributions |  [link](https://github.com/mitmath/18338/blob/master/projects/2023/bryutkinandrey_LATE_152454_3887228_Painleve%CC%81%20Systems%20Presentation-1.pdf) |  [link](https://github.com/mitmath/18338/blob/master/projects/2023/bryutkinandrey_LATE_152454_3887229_random_matrices_final_project-2-1.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/bryutkinandrey_LATE_152454_3887230_testing_painleve-2.m)|
| |Bünyamin Kartal & Maison Clouâtré| Random Quantum Density Operators | | [link](https://github.com/mitmath/18338/blob/master/projects/2023/clouatremaison_118505_3875544_RandomQuantumDensityOperators%E2%80%94Final.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/clouatremaison_118505_3875545_Code1.ipynb) [link](https://github.com/mitmath/18338/blob/master/projects/2023/clouatremaison_118505_3875546_Code2.ipynb)|
| |Matthew Esmaili Mallory| Free Probability & The Free Central Limit Theorem | [link](https://github.com/mitmath/18338/blob/master/projects/2023/esmailimatthew_LATE_163052_3886901_Final%20Prez%20338.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/esmailimatthew_163052_3873370_Final_Project___18_338%20(1).pdf)| |
| |Mali Halac| Determinantal Point Processes (DPPs) Improve KNN Clasifier Performance in Bioimpendence Analysis | | [link](https://github.com/mitmath/18338/blob/master/projects/2023/halacmali_163226_3875472_MIT_18_338_Final_Paper.pdf) | |
| |Tanshiq Kumar| Random Matrix Theory and Generalization in Neural Networks | [link](https://github.com/mitmath/18338/blob/master/projects/2023/kumartanishq_114210_3876069_338_Final_Project_Slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/kumartanishq_114210_3876068_TK_338_Final_Report_Complete.pdf)| |
| |Yizhou Liu| Are Interactions Real | | [link](https://github.com/mitmath/18338/blob/master/projects/2023/liuyizhou_134814_3871872_18_338_project_final.pdf)| |
| |Shyam Narayan| Tail Bounds on the Smallest Singular Value of a Rectangular Random Matrix| [link](https://github.com/mitmath/18338/blob/master/projects/2023/NarayananShyamSlides.pptx) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/narayananshyam_86941_3874130_18.338%20Final%20Project.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2023/narayananshyam_86941_3874131_18.338%20Final%20Project%20Experiments.html) |
| |Nicholas Stiles| Eigenvectors of the Correlation Matrix | [link](https://github.com/mitmath/18338/blob/master/projects/2023/stilesnicholas18_338_project_presentation.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2023/stilesnicholas_LATE_120061_3879963_18.338%20project%20report.pdf)| |
| |Songchen Tan| Application of RMT in some Discrete Optimization Problems | | [link](https://github.com/mitmath/18338/blob/master/projects/2023/tansongchen_97375_3874845_Application%20of%20RMT%20in%20some%20Discrete%20Optimization%20Problems.pdf) | |
| | Harry Walden | Roots of Random Polynomials with Integral Geometry | [link](https://github.com/mitmath/18338/blob/master/projects/2023/Roots%20of%20random%20Polynomials%20with%20integral%20geometry.pptx) | [link](https://github.com/mitmath/18338/blob/master/projects/2023/waldenharry_154930_3866094_hjw57_18338_report.pdf) | [link](https://github.com/mitmath/18338/tree/master/projects/2023/code)|
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
| 2017 | Sungwoo Jeong| Jacobians of Matrix decompositions | | [link](https://github.com/mitmath/18338/blob/master/projects/2017/JacobianofMatrixDecompositions.pdf) | |
| 2016 |Antoni Musolas |Differential geometrical approach to covariance estimation | | | |
| |Manishika Agaskar |	Diffraction gain of free-space optical communications in atmospheric turbulence | | | |
| |Hong Hu | Spectral Initialization and its Performance Analysis| | | |
| |Anuran Makur	|Maximal Correlation Functions: Hermite, Laguerre, and Jacobi  | | | |
| | John Urschel |On the Minimal Eigenpair of Erdos-Renyi Graphs | | | |
| |John Napp |RMT and the complexity of linear optics | | | |
| |Alex Wein |Random Matrix Contiguity | | | |
| |Ravi Bajaj |Central Limit Theorem for log-Determinant of Wigner Matrices | | | |
| |Brandon Tran |Application of phase transitions in spiked covariance matrices to synchronization problems | | | |
|2012|Chenhui Hu | Spectral Perturbation of Small-World Networks | [link](https://github.com/mitmath/18338/blob/master/projects/2012/ch_slides.pdf)| [link](https://github.com/mitmath/18338/blob/master/projects/2012/ch_report.pdf) |  |
| | Mina Karzand | Random Matrix Theory and Non-coherent MIMO Communications | [link](https://github.com/mitmath/18338/blob/master/projects/2012/mk_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/mk_report.pdf) | | 
| |Charlotte Kiang | 	RMT Applications to Control Theory | [link](https://github.com/mitmath/18338/blob/master/projects/2012/ck_slides.pptx) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/ck_report.pdf) | |
| |Yi Sun| Virasoro Constraints and the Tracy-Widom Law | [link](https://github.com/mitmath/18338/blob/master/projects/2012/ys_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/ys_report.pdf) | | 
| |Matt Welborn | Density of states of disordered systems via free addition  | [link](https://github.com/mitmath/18338/blob/master/projects/2012/mw_slides.pptx) | | | 
| | Helen Xie | RMT Applications to Transport Property | [link](https://github.com/mitmath/18338/blob/master/projects/2012/hx_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/hx_report.pdf) | | 
| | Yufei Zhao| Spectral Distributions of Random Graphs | [link](https://github.com/mitmath/18338/blob/master/projects/2012/yz_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/yz_report.pdf)| |
| |Yi Zeng |North Pole Problem | [link](https://github.com/mitmath/18338/blob/master/projects/2012/zy_slides.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2012/zy_report.pdf) | |
| 2009| Gregory Minton | Moments of Random Orthogonal Matrices | [link](https://github.com/mitmath/18338/blob/master/projects/2009/final-report.pdf) | [link](https://github.com/mitmath/18338/blob/master/projects/2009/momentinterface.html)|
 
