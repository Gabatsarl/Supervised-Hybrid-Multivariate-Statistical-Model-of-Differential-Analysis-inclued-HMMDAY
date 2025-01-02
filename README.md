# Introduction

L'analyse différentelle (DA) permet d'étudier les gènes différentielles en fonction des conditions d'expérience. Si nous notons $X_g$ l'abondance des gènes et C les conditions métagénomiques, le DA s'exprime par

$$
X_g=C \beta+\varepsilon \cdot
$$

Pour chaque gène $g$, nous étudions sa significativité par la valeur de la pvaleur associée au paramètre $\beta$. Cependant, lorsque nous voulons faire cette analyse différence en incluant une variable d'intérêt $Y$, la problématique devient différent et implique donc 

$$Y=X_g \psi+C \beta+\varepsilon \cdot$$

Cette nouvelle problématique devient relativement compte tenu du faite que c'est une régression mais des difficultés comme la taille des échantillons métagénomique, de la sparsité des données métagénomiques limites l'application des méthodes de régression classique. 

# Méthode proposée
Nous proposons l'utilisation d'une stratégie classique issu de la méthode de l'équation $Y=X_g \psi+C \beta+\varepsilon$. Nous avons calculé les metrics comme le $\lam et tracé les graphiques de qqplot pour qualifier la qualité de l'ajustement avant de conclure.


# Écriture du modèle 
Le modèle de l'équation proposée est $Y=X_g \psi+C \beta+\varepsilon$ permet d'analyser l'éfficacité de chaque gène $g$ par rapport à la condition $C$ sur la variable $Y$. 

# Études de simulation
Dans cette étude de simulation, nous simulons deux types de données $X$ : les matrice $X$ est générée soit à partir des données observées $X_{obs}$ ou soit de manière aléatoire.

1. Simulation à partir des données existantes. Cette simulation prend en "entrées" des paramètres estimées à partir des données observ.es . 

2. Simulation sans les données. Cette simulation prend en "entrées" des paramètres aléatoirement choisis par l'utilisateur pour générer la matrice de données $X$.

# Résultats

Le modèle comparative  a permis de détecter le modèle qui s'ajuste le mieux et d'obtenir les résultats souhaités.
Cependant, une amélioration et une extension de ce modèle en multidimension est souhaitable.



