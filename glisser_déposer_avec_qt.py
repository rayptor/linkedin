import sys
import uuid
from PyQt6.QtCore import Qt, QIODevice, QDataStream
from PyQt6.QtWidgets import (
    QWidget, QApplication, QGridLayout, QVBoxLayout, QHBoxLayout,
    QLabel, QListWidget, QListWidgetItem, QPushButton, QFrame
)

PRIX_PAR_DEFAUT = {
    "PIZZA": {
        "Reine": 12.0, "Napolitaine": 13.0, "Marguerite": 14.0,
        "Cannibale": 19.0, "Norvégienne": 22.0, "Truffe": 25.0
    },
    "BOISSONS": {"Eau plate": 1.5, "Eau gazeuse": 2.0, "Jus de fruits maison": 4.5},
    "DESSERTS": {"Tarte aux pommes": 4.5, "Muffin au chocolat": 4.0, "Tropézienne": 5.5, "Tiramisu": 6.5}
}

class ListeArticles(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(False)
        self.setDragDropMode(QListWidget.DragDropMode.DragOnly)

    def dragEnterEvent(self, event):
        event.ignore()

    def dropEvent(self, event):
        event.ignore()


class ListeCommandeArticles(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.setDragDropMode(QListWidget.DragDropMode.DropOnly)
        self.itemDropped = None
        self.parent_app = parent

        self.setStyleSheet("""
            QListWidget {
                background-color: #FFFFFF;
                border: 2px solid #FCC;
                border-radius: 10px;
                padding: 5px;
            }
        """)

    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.setDropAction(Qt.DropAction.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        mime_data = event.mimeData()
        if mime_data.hasFormat("application/x-qabstractitemmodeldatalist"):
            encoded_data = mime_data.data("application/x-qabstractitemmodeldatalist")
            data_stream = QDataStream(encoded_data, QIODevice.OpenModeFlag.ReadOnly)
            data_stream.setVersion(QDataStream.Version.Qt_6_0)

            while not data_stream.atEnd():
                _ = data_stream.readInt()
                _ = data_stream.readInt()
                map_items = data_stream.readInt()

                for _ in range(map_items):
                    cle = data_stream.readInt()
                    valeur = data_stream.readQVariant()

                    if Qt.ItemDataRole(cle) == Qt.ItemDataRole.UserRole:
                        categorie, libelle, prix = valeur
                        prix_article = PRIX_PAR_DEFAUT.get(categorie, {}).get(libelle, prix)
                        texte_formate = f"{libelle} - {prix_article:.2f} €"
                        article_depose = QListWidgetItem(texte_formate)
                        article_a_commander = (categorie, libelle, prix_article)
                        if categorie == "PIZZA":
                            pizza_identifiant = str(uuid.uuid4()) # ID aléatoire
                            article_a_commander = (categorie, libelle, prix_article, pizza_identifiant)
                            self.parent_app.dernier_element_ajoute = article_depose
                        article_depose.setData(Qt.ItemDataRole.UserRole, article_a_commander)
                        self.addItem(article_depose)
                        if self.itemDropped:
                            self.itemDropped()
            event.acceptProposedAction()
        else:
            event.ignore()

class PizzeriaBase(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Commande Pizzeria : démo glisser-déposer")
        self.prix = PRIX_PAR_DEFAUT
        self.initialiser_interface()
        self.setStyleSheet(self.style_global())

    def style_global(self):
        return """
            QWidget {
                background-color: #FF6347;
                color: #333;
                font-family: Arial, sans-serif;
            }
            QLabel {
                color: #FFFFFF;
            }
            QLabel#montant_libelle {
                font-size: 20px;
                font-weight: bold;
                color: #FFFFFF;
            }
            QListWidget {
                background-color: #FFFFFF;
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 5px;
                color: #333;
            }
            QListWidget:hover {
                background-color: #FF8569;
                border: 1px solid #fff;
                border-radius: 5px;
            }
            QListWidget::item {
                padding: 3px;
            }
            QListWidget::item:selected {
                background-color: #FF6347;
                color: white;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                padding: 10px 20px;
                border: none;
                border-radius: 5px;
                font-size: 18px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
                border: 1px solid #fff;
                border-radius: 5px;
            }
            QFrame {
                background-color: transparent;
            }
            QFrame[frameShape="4"] {
                background-color: #FFDAB9;
                height: 2px;
            }
        """

    def creer_liste_par_categorie(self, titre, articles):
        layout = QVBoxLayout()
        layout.addWidget(QLabel(f"<b>{titre}</b>"))
        liste_widgets_articles = ListeArticles()
        liste_widgets_articles.setDragEnabled(True)

        for nom in articles:
            prix = PRIX_PAR_DEFAUT[titre][nom]
            article = QListWidgetItem(f"{nom} - {prix:.2f} €")
            article.setData(Qt.ItemDataRole.UserRole, (titre, nom, prix))
            article.setFlags(article.flags() | Qt.ItemFlag.ItemIsDragEnabled)
            liste_widgets_articles.addItem(article)
        layout.addWidget(liste_widgets_articles)
        return layout, liste_widgets_articles

    def initialiser_interface(self):
        layout_base = QHBoxLayout()
        self.listes = {}
        categorie_prix = list(PRIX_PAR_DEFAUT.keys()) 
        
        layout_categorie = QVBoxLayout()
        layout_grille = QGridLayout()
        ligne = 0
        colonne = 0
        for categorie in categorie_prix:
            layout, liste_widgets = self.creer_liste_par_categorie(categorie, self.prix[categorie])
            layout_grille.addLayout(layout, ligne, colonne)
            self.listes[categorie] = liste_widgets
            colonne += 1
            if colonne > 2:
                colonne = 0
                ligne += 1
        layout_categorie.addLayout(layout_grille)
        layout_base.addLayout(layout_categorie)

        layout_commande = QVBoxLayout()
        layout_commande.addWidget(QLabel("<b>Commande</b>"))

        self.liste_commande = ListeCommandeArticles(self)
        self.liste_commande.itemDropped = self.mettre_a_jour_montant
        self.liste_commande.itemClicked.connect(self.dernier_element_clique)
        layout_commande.addWidget(self.liste_commande)

        self.montant_libelle = QLabel("Montant : 0.00 €")
        self.montant_libelle.setObjectName("montant_libelle")
        self.montant_libelle.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        layout_commande.addWidget(self.montant_libelle)

        separateur = QFrame()
        separateur.setFrameShape(QFrame.Shape.HLine)
        separateur.setFrameShadow(QFrame.Shadow.Sunken)
        layout_commande.addWidget(separateur)

        bouton_annuler = QPushButton("Annuler la commande")
        bouton_annuler.clicked.connect(self.annuler_commande)
        layout_commande.addWidget(bouton_annuler)
        layout_commande.addStretch(1)
        layout_base.addLayout(layout_commande)
        self.setLayout(layout_base)
        self.setFixedSize(850, 400)

    def dernier_element_clique(self, item):
        item_data = item.data(Qt.ItemDataRole.UserRole)
        if item_data and item_data[0] == "PIZZA":
            self.dernier_element_ajoute = item

    def annuler_commande(self):
        self.liste_commande.clear()
        self.dernier_element_ajoute = None
        self.mettre_a_jour_montant()

    def mettre_a_jour_montant(self):
        total = 0.0
        for i in range(self.liste_commande.count()):
            article = self.liste_commande.item(i)
            article_data = article.data(Qt.ItemDataRole.UserRole)
            if article_data and isinstance(article_data, tuple) and len(article_data) >= 2:
                category = article_data[0]
                name = article_data[1]
                prix_article_ajoute = PRIX_PAR_DEFAUT.get(category, {}).get(name, 0.0)
                total += prix_article_ajoute

        self.montant_libelle.setText(f"Montant : {total:.2f} €")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    fenetre = PizzeriaBase()
    fenetre.show()
    sys.exit(app.exec())