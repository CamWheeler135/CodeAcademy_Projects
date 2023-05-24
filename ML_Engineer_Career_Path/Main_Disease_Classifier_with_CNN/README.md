# Can we use Your Immune Data to Predict Cancer?

### ❗️Project Currently Being Built❗️

In this project, I aim to use immune cell data in order to predict the variant of cancer a patient has. I am going to be quite informal throughout the project and assume little to no previous knowledge in either medical biology or machine learning on behalf of the reader. Any technical topics will be explained so anybody can pick up this README and understand what I am doing, how I am doing it and why! If after reading this you want to delve deeper into the subject, I have linked plenty of papers in [this section](#papers-and-interesting-articles). 

## TLDR ⚠️

Being able to accurately diagnose cancer is crucial for appropriate and effective treatment. Current methods used in cancer diagnosis look at imaging, biopsies and lab tests. In the age of next generation sequencing (NGS) we can quickly and cost effectively sequence cells producing HUGE amounts of data that we can analyze using machine learning methods. Our adaptive immune system is what allows us to combat disease that our innate immune system would not recognize 🦠, a key part of this is the T cell. This project aims to use a deep learning model (CNN) to see if we can classify the specific cancer a patient has presented with based on their repertoire data. Being able to do such a task could lead to another tool in the clinicians arsenal when it comes to diagnosing patients. 

## Jump Right In ⤵️

- [Introduction](#introduction)
- [Aims](#aims)
- [Data](#data)
- [Findings](#findings)
- [Methods](#methods)
- [Papers and Interesting Articles](#papers-and-interesting-articles)

## Introduction
#### *More than the TLDR.* 🔎


## Aims
#### *What are my aims for the project?* 🎯
- Understand the workings of cancer diagnosis and the pros and cons of current methods. 
- Further investigate immune data and its usefulness in building clinical machine learning models.
    - Increase my ability to visualize and communicate findings with regards to medical data and machine learning. 
- Increase my understanding of end to end project building. 
- Learn PyTorch fundamentals. 
- Learn to communicate the technical methods and knowledge properly to individuals outside of the field.

## Data
#### *What data is being used in the project?* 📂

The data being used in the project is coming from the TCR data base [(TCRdb)](http://bioinfo.life.hust.edu.cn/TCRdb/#/browse).

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

It is worth noting that there is some class imbalance here, some cancers contain a larger number of repertoires than others so we will have to deal with this in the preprocessing stage to ensure that our model is effective. 

## Findings
#### *What did we discover and learn in the project?* 🔬
- This is where our findings will go. 

## Methods
#### *What did we use in the project?* 🧪
- This is where our methods will go.

## Papers and Interesting Articles
#### *Formal reading and interesting articles I used to build the project.* 📚
- This is where our references. 
