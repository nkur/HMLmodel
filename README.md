---
<div align='center'>
# Human Motor Learning Dynamics in High-dimensional Tasks
Ankur Kamboj, Rajiv Ranganathan, Xiaobo Tan, Vaibhav Srivastava | 2024  
[![Paper](https://img.shields.io/badge/PLOS_CB-2024-red)](https://doi.org/10.1371/journal.pcbi.1012455)

</div>

_TL;DR_: This work models human skill learning in high-dimensional motor spaces leveraging the concept of motor synergies, internal model theory of motor learning, and adaptive control.  
## Summary
Conventional approaches to enhance movement coordination, such as providing instructions and visual feedback, are often inadequate in complex motor tasks with multiple degrees of freedom (DoFs). To effectively address coordination deficits in such complex motor systems, it becomes imperative to develop interventions grounded in a model of human motor learning; however, modeling such learning processes is challenging due to the large DoFs. In this paper, we present a computational motor learning model that leverages the concept of motor synergies to extract low-dimensional learning representations in the high-dimensional motor space and the internal model theory of motor control to capture both fast and slow motor learning processes. We establish the model’s convergence properties and validate it using data from a target capture game played by human participants. We study the influence of model parameters on several motor learning trade-offs such as speed-accuracy, exploration-exploitation, satisficing, and flexibility-performance, and show that the human motor learning system tunes these parameters to optimize learning and various output performance metrics.  

---

This repository contains MATLAB implementation of the Human Motor Learning model presented in the corresponding paper.

Additional files are required to reproduce the results presented in the paper concerning various human motor learning behavior.

The folder `\pierellaModelFits` contains MATLAB codes that compare the HML model with the state-of-the-art human motor learning dynamics model presented in [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1007118).
