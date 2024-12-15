# Introduction

L'analyse différentelle (DA) permet d'étudier les gènes différentielles en fonction des conditions d'expérience. Si nous notons $X_g$ l'abondance des gènes et C les conditions métagénomiques, elle s'exprime par

$$
X_g=C \beta+\varepsilon \cdot
$$

Pour chaque gène $g$, nous étudions sa significativité par la valeur de la pvaleur associée au paramètre $\beta$. Cependant, lorsque nous voulons faire cette analyse différence en incluant une variable d'intérêt $Y$, la problématique devient différent et implique donc 

$$Y=X_g \psi+C \beta+\varepsilon \cdot$$

Cette nouvelle problématique devient relativement compte tenu du faite que c'est une régression mais des difficultés comme la taille des échantillons métagénomique, de la sparsité des données métagénomiques limites l'application des méthodes de régression classique. Nous proposons l'utilisation d'une stratégie classique issu de la méthode de l'équation \eqre{equation1}.
