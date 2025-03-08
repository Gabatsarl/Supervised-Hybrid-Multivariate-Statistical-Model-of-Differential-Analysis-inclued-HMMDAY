# Introduction

L'analyse différentielle (DA) permet d'étudier les gènes différentiellement exprimés en fonction des conditions expérimentales. Si nous notons $X_g$, l'abondance des gènes et $C$ les conditions métagénomiques, le DA s'exprime par :

$$
X_g=C \beta+\varepsilon,
$$

where $\beta$ représente les effets fixes et $\varepsilon$, les erreurs aléatoires.

Pour chaque gène $g$, nous étudions sa significativité par la pvaleur associée au paramètre $\beta$. 

Cependant, lorsque nous voulons faire cette analyse différentielle en incluant une variable d'intérêt $Y$, la problématique devient différente et implique donc l'utilisation de l'équation :

$$Y=X_g \psi+C \beta+\varepsilon,$$

Où :
- $Y$ représente la variable dépendante (ou la matrice de données métagénomiques observées),
- $X_g$ est l'abondance du gène $g$,
- $\psi$ et $\beta$ sont des coefficients,
- $C$ est une condition binaire (par exemple, une condition expérimentale),
- $\varepsilon$ est l'erreur du modèle, supposée suivre une distribution **normale** ou de **Student** selon le cas.

Cette équation est relativement simple compte tenu du faite que c'est une régression linéaire. Malheureusement des difficultés comme la taille des échantillons métagénomique, de la sparsité des données métagénomiques, etc., limitent l'application des méthodes de régression classique. 

# Méthode proposée
Nous proposons l'utilisation d'une stratégie classique issu de la méthode de l'équation $Y=X_g \psi+ C \beta+\varepsilon$. Nous avons calculé les metrics comme le $\lambda$ et tracé les graphiques de QQplot pour qualifier l'ajustement avant de conclure. Cet ajustement de modèle à partir de DA incluant $Y$, nous la nommons le modèle ***HMMDAY***.

# Écriture du modèle
Le modèle de l'équation proposée est $Y=X_g \psi+C \beta+\varepsilon$ permet d'analyser l'éfficacité de chaque gène $g$ par rapport à la condition $C$ sur la variable $Y$. 

# Études de simulation du modèle ***HMMDAY***

## Contexte

Dans le cadre de cette étude de simulation, l'objectif principal est de générer des matrices de données $X$ et des erreurs $\epsilon$ sous différentes conditions. Nous simulons la matrice de données $X$ selon deux stratégies distinctes : 

1. **Simulation à partir des données observées**  $X_{\text{obs}}$
2. **Simulation aléatoire de $X$**

En parallèle, nous faisons l'hypothèse que les erreurs suivent une distribution normale dans un cas et une distribution de Student dans un autre. Enfin, nous supposons que la condition $C$ est binaire.

## Objectifs globaux

L'objectif de cette étude de simulation est de :
- Tester l'applicabilité du modèle **HMMDAY** dans des contextes réels et généralisés.
- Analyser comment différentes distributions d'erreurs et types de données affectent les performances du modèle.
  
Nous espérons que les résultats obtenus grâce à cette approche de simulation aideront à affiner et à adapter le modèle ***HMMDAY*** aux défis complexes des données métagénomiques.

## Stratégies de simulation

### 1. Simulation de $X$ à partir des données observées $X_{\text{obs}}$

Cette stratégie utilise les données observées $X_{\text{obs}}$ comme base pour générer la matrice de données $X$. Le processus de simulation repose sur les **paramètres estimés** à partir des données observées. Cette approche permet de reproduire les variabilités biologiques, techniques et compositionnelles caractéristiques des données métagénomiques existantes. Elle est particulièrement utile pour étudier l'applicabilité du modèle sur des **données réelles**.

**Objectifs de cette simulation :**
- Reproduire les caractéristiques des données réelles.
- Étudier comment le modèle ***HMMDAY*** réagit face à des données observées, ce qui est essentiel pour valider son efficacité dans des situations réelles.

### 2. Simulation de $X$ sans données existantes $X_{\text{obs}}$

Dans cette approche, la matrice de données $X$ est générée à partir de **paramètres aléatoires**, choisis par l'utilisateur. Cette stratégie est utilisée pour simuler des ensembles de données d’ordre général, sans lien direct avec des données observées. Cela permet de tester l’applicabilité du modèle ***HMMDAY*** sur des données **synthétiques**, et ainsi d’étudier sa robustesse et sa capacité d’adaptation à des situations plus générales ou inconnues.

**Objectifs de cette simulation :**
- Tester le modèle ***HMMDAY*** sur des données générées aléatoirement.
- Évaluer la performance du modèle dans des contextes variés où les données réelles ne sont pas disponibles.

## Hypothèses

1. Les erreurs $\epsilon$ suivent une distribution **normale** dans le premier cas.
2. Les erreurs $\epsilon$ suivent une distribution de **Student** dans le deuxième cas.
3. La condition $C$ est supposée être **binaire**.

Ces hypothèses sont importantes pour la validation des résultats des simulations et permettent d’étudier le comportement du modèle ***HMMDAY*** dans différents environnements statistiques.

# Résultats

Le modèle comparative  a permis de détecter le modèle qui s'ajuste le mieux et d'obtenir les résultats souhaités. Cependant, une amélioration et une extension de ce modèle en multidimension est souhaitable.


# Structure du programmes

## Programme "clean.R"
Ce programme permet de nettoyer les données d'**abondance métagénomique** et d'**annotation** disponibles. Nous disposons de plusieurs jeux de données.

## programme "umap.R"
Pour ce programme, nous mettons en oeuvre la méthode de réduction de dimension pour obtenir les deux vecteurs des axes (Axe 1 et Axe 2). L'information obtenue sera utiliser dans les régressions individuelles. Le fichier **umap.R** utile les résultats du fichier **clean.R**.

## Programme "slide.R"
Pour ce programmes, nous décomposons la matrice **X** en 100 jeux de données ayant le **n** lignes et p/100 colonnes. L'objectif de cette subdivision est de maximiser les ressources de **Compute Canada**. Le fichier **slide.R** utile les résultats du fichier **clean.R**

## Programme "method.R"
Pour ce programme, nous mettons en oeuvre le modèle **HMMDAY** sur le jeu de données de **slide.R** 

## Programme 
