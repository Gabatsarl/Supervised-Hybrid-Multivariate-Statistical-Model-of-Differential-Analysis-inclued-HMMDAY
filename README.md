# Introduction

L'analyse différentielle (DA) permet d'étudier les gènes différentielle exprimés en fonction des conditions expérimentales. Si nous notons $X_g$ l'abondance des gènes et $C$ les conditions métagénomiques, le DA s'exprime par :

$$
X_g=C \beta+\varepsilon,
$$

where $\beta$ représente les effets fixes et $\varepsilon$, les erreurs aléatoires.

Pour chaque gène $g$, nous étudions sa significativité par la pvaleur associée au paramètre $\beta$. 

Cependant, lorsque nous voulons faire cette analyse différentielle en incluant une variable d'intérêt $Y$, la problématique devient différente et implique donc l'utilisation de l'équation :

$$Y=X_g \psi+C \beta+\varepsilon \cdot$$

Cette équation est relativement simple compte tenu du faite que c'est une régression linéaire. Malheureusement des difficultés comme la taille des échantillons métagénomique, de la sparsité des données métagénomiques, etc limitent l'application des méthodes de régression classique. 

# Méthode proposée
Nous proposons l'utilisation d'une stratégie classique issu de la méthode de l'équation $Y=X_g \psi+ C \beta+\varepsilon$. Nous avons calculé les metrics comme le $\lambda$ et tracé les graphiques de QQplot pour qualifier la qualité de l'ajustement avant de conclure. Cet ajustement de modèle à partir de DA incluant $Y$, nous la nommons le modèle $\textbf{HMMDAY}$.


# Écriture du modèle
Le modèle de l'équation proposée est $Y=X_g \psi+C \beta+\varepsilon$ permet d'analyser l'éfficacité de chaque gène $g$ par rapport à la condition $C$ sur la variable $Y$. 

# Études de simulation
Dans cette étude de simulation, nous avons besoin de la matrice de données $X$ et les erreurs $\varepsilon$. Pour ce fait, nous simulons la matrice de données $X$ par deux stratégies : les matrice $X$ est générée soit à partir des données observées $X_{obs}$ ou soit de manière aléatoire. Nous supposons que les erreurs suivent une distribution normale d'une part et une distribution de Student d'autre part. Enfin, nous suppons que la condition $C$ est binaire.

1. Simulation de $X$ à partir des données existantes $X_{obs}$. Cette simulation prend en "entrées" des paramètres estimées à partir des données observées. Cette simulation permet d'avoir une replication des variabilités biologiques, technique et compositionnele des données métagénomique. Ceci permet d'étudier l'applicabilité du modèle sur des données réelles.

2. Simulation de $X$ sans les données existantes $X_{obs}$. Cette simulation prend en "entrées" des paramètres aléatoirement choisis par l'utilisateur pour générer la matrice de données $X$. Pour étudier, l'applicabilité du modèle HMMDAY sur des données d'ordre général 

# Résultats

Le modèle comparative  a permis de détecter le modèle qui s'ajuste le mieux et d'obtenir les résultats souhaités.
Cependant, une amélioration et une extension de ce modèle en multidimension est souhaitable.



