# Can we use Your Immune Data to Predict Cancer?

### ❗️Project Currently Being Built❗️

In this project, I aim to use immune cell data in order to predict the variant of cancer a patient has. This in an informal summary of what the project is and what happened throughout. I assume little to no previous knowledge in either medical biology, programming or machine learning on behalf of the reader. Any technical topics will be explained so anybody can pick up this README and understand what I am doing, how I am doing it and why! After reading this page you may want to delve deeper into the subject, I have linked plenty of papers in [this section](#papers-and-interesting-articles). 

## TLDR ⚠️

Being able to accurately diagnose cancer is crucial for appropriate and effective treatment. Current methods used in cancer diagnosis look at imaging, biopsies and lab tests. In the age of next generation sequencing (NGS) we can quickly and cost effectively sequence cells producing HUGE amounts of data that we can analyze using machine learning methods. Our adaptive immune system is what allows us to combat disease that our innate immune system would not recognize 🦠, a key part of this is the T cell. This project aims to use a deep learning model (CNN) to see if we can classify the specific cancer a patient has presented with based on their repertoire data. Being able to do such a task could lead to another tool in the clinicians arsenal when it comes to diagnosing patients. 

## Jump Right In ⤵️

- [Introduction](#introduction)
    - [Cancer Biology](#cancer-introduction)
    - [The Immune System](#introduction-into-the-immune-system)
    - [Linking Cancer and the CDR3$\beta$ Sequence](#linking-cancer-and-cdr3-sequence)
    - [Machine Learning](#introduction-into-machine-learning)
- [Aims](#aims)
- [Data](#data)
- [Findings](#findings)
- [Methods](#methods)
- [Papers and Interesting Articles](#papers-and-interesting-articles)

## Introduction
#### *More than the TLDR.* 🔎

#### Cancer Introduction. 
- What is it?
- How do we detect it?
- How do we treat it?

#### Introduction into the immune system. 

As we go about our daily lives, our bodies are constantly coming into contact with immunogens (things that illicit an immune response). This could come from your colleague sneezing in a meeting, taking a walk in the park with hay fever or like me, you chew your nails; the air we breath and the surfaces we touch are FULL of bad things. I know right, the world is kinda nasty 🦠. Luckily for us, our bodies have lots of defenses to make sure we stay alive without having to shower in Dettol every time we get back from work or school.

This is where our *immune system* comes in. There are two main systems you need to understand, the *innate* system and the *adaptive* system. The innate immune system contains all the physical and chemical barriers like our skin, mucosa tract and digestive system. These take the approach of, "If it can't get in, it can't hurt us", its rather blunt, yet still effective. There also exists cell mediated innate immunity which include cells like phagocytes that engulf and digest invaders, presenting the remains on their cell surface (the fact they can do this is REALLY important so remember this for later). The problem is that these innate immune cells cannot recognize all the invaders that enter our bodies. Our *adaptive* immune system covers this; it is able to specifically target a particular immunogen that is causing the problem, while in the processes *remember* what that immunogen looked like, just in case it comes back later (this 'memory' mechanism is what immunologists aim to exploit with vaccines). The adaptive immune system is made by of T lymphocytes and B lymphocytes, but we are going to focus on the T lymphocyte (T cell). A T cell can have several types of jobs depending on its variant, it can either go out and fight the infection, help recruit other immune cells to the fight or memorize the immunogen for later. 

Now that we are familiar with the role of the adaptive immune system and the T cell, we need to talk about how to T cell recognizes these invaders. Each T cell has its own T cell receptor (TCR), the majority of TCRs in humans consist of $\alpha$ and $\beta$ subunits. On the $\beta$ subunit exists a specific sequence called the Complementary-Determining-Region-3-Beta (CDR3β), this sequence is what is hypothesized by immunologists to control what the T cell can bind to. If a TCR is thought to bind to a particular immunogen, it is said to have a high *affinity* to that particular sequence. When T cells are being developed, the TCR CDR3$\beta$ sequence is formed through V(D)J recombination. This is a complex process that takes the T cell DNA and chops it up in random places, this process (along with some errors that are usually made when stitching the sequence back together) is though to produce over $10^{15}$ possible CDR3$\beta$ sequences, that's 1,000,000,000,000,000 possible sequences❗️ The fact that we can do this, allows us to defend ourselves so well against so many pathogens that we come across in our daily activities. 

#### Linking Cancer and CDR3β Sequence

- Introduction into the machine learning. 
    - What are we using?
    - What task are we trying to do?

#### Introduction into Machine Learning


## Aims
#### *What are my aims for the project?* 🎯

- Understand the workings of cancer diagnosis and the pros and cons of current methods. 
- Further investigate immune data and its usefulness in building clinical machine learning models.
    - Increase my ability to visualize and communicate findings with regards to medical data and machine learning. 
     Learn to communicate the technical methods and knowledge properly to individuals outside of the field.
- Increase my understanding of end to end project building. 
- Learn PyTorch fundamentals. 

## Data
#### *What data is being used in the project?* 📂

The data being used in the project is coming from the TCR data base [(TCRdb)](http://bioinfo.life.hust.edu.cn/TCRdb/#/browse). It is build of TCR CDR3β sequences taken from a several solid tumour cancers and blood cancers in humans. This data also comes along with V, D, and J region information that we are going to use in the exploratory data analysis. Importantly, the database has already done a lot of the cleaning for us, but I have build a data cleaning module anyway just in case. It will also be useful when trying to classify samples that have not come from the database as we can handle and process them properly. These samples are taken from a range of experiments, patients and cell sources, furthermore, the times in which these samples were taken also differ. Although this might not seem ideal as there doesn't appear to be much homogeny in the data collection methods. THIS IS GREAT! Medical data is never perfect, we can never get samples at the same time, from patients who are the same age, with the same disease progress, from the same source. Having a model that can perform well with diverse data has a much better chance in the clinic than one that is really sensitive to such issues. 

The cancers that are included in version 1 of the project are:

- Esophageal
- Breast
- Gastric
- Liver
- Lung 
- Glioblastoma 
- Ovarian
- Head and Neck
- Colorectal
- Pancreatic

## Methods
#### *What did we use in the project?* 🧪

#### Data Cleaning
- What we are doing to clean the data. 

## Findings
#### *What did we discover and learn in the project?* 🔬
- This is where our findings will go. 

## Papers and Interesting Articles
#### *Formal reading and interesting articles I used to build the project.* 📚
- This is where our references. 
