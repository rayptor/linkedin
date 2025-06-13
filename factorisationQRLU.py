# fichier factorisationQRLU.py

import sys
import numpy as np
import scipy.linalg as spla
from PyQt6.QtCore import Qt, QLocale, QCoreApplication, QRegularExpression
from PyQt6.QtGui import QFont, QRegularExpressionValidator
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSpinBox,
    QRadioButton, QLabel, QPushButton, QLineEdit, QGridLayout, QMessageBox,
    QGroupBox, QSizePolicy, QButtonGroup, QFrame, QFileDialog
)

class FactorisationQRLU(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.ordre: int = 3
        self.decomposition: str = "QR"
        self.entrees: list = []
        self.libelles_q: list = []
        self.libelles_r: list = []
        self.libelles_l: list = []
        self.libelles_u: list = []

        self.setWindowTitle("Factorisation QR / LU")
        self.setFixedSize(750, 800)
        self.setLayout(QVBoxLayout())
        self.setStyleSheet("""
           QWidget {
                background-color: #2E2E2E;
                color: #E0E0E0;
                font-size: 14px;
            }
            QPushButton {
                background-color: #3A3A3A;
                color: white;
                border: 1px solid #4F4F4F;
                border-radius: 6px;
                padding: 10px 20px;
                font-size: 16px;
                font-weight: bold;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #4A4A4A;
                border-color: #5F5F5F;
            }
            QPushButton:pressed {
                background-color: #2F2F2F;
                border-color: #3A3A3A;
            }
            QPushButton:focus {
                border: 2px solid #00BFFF;
            }
            QLabel {
                font-size: 15px;
                color: #C0C0C0;
                font-weight: bold;
            }
            QRadioButton {
                background-color: #383838;
                color: #E0E0E0;
                font-size: 15px;
                font-weight: normal;
                padding: 5px 0;
                margin-right: 15px;
            }
            QRadioButton::indicator {
                width: 18px;
                height: 18px;
                border-radius: 9px;
                border: 2px solid #606060;
                background-color: #4A4A4A;
                margin-right: 8px;
            }
            QRadioButton::indicator:hover {
                border: 2px solid #00BFFF;
            }
            QRadioButton::indicator:checked {
                background-color: #00BFFF;
                border: 2px solid #00BFFF;
            }
            QRadioButton::indicator:checked:hover {
                background-color: #00A3D9;
            }
            QLineEdit {
                background-color: #3A3A3A;
                color: #F0F0F0;
                font-size: 15px;
                border-radius: 4px;
                border: 1px solid #555555;
                padding: 2px;
                text-align: center;
            }
            QLineEdit:hover {
                border: 1px solid #777777;
            }
            QLineEdit:focus {
                border: 2px solid #00BFFF;
            }
            QSpinbox {
                font-size: 15px;
                color: #C0C0C0;
                font-weight: bold;
            }
            QSpinBox:hover
            {
                border: 1px solid #00BFFF;
            }
            QSpinBox::down-button {
                subcontrol-position: left;
                height: 16px;
                margin: 0 5px 1px 0;
                padding: 0 1px 1px 1px;
            }
            QSpinBox::up-button {
                subcontrol-position: right;
                height: 16px;
                margin: 0 5px 1px 0;
                padding: 0 1px 1px 1px;
            }
            QGroupBox {
                background-color: #383838;
                border: 1px solid #505050;
                border-radius: 10px;
                margin-top: 25px;
                font-size: 16px;
                font-weight: bold;
                color: yellow;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top center;
                padding: 0 10px;
                background-color: #303030;
                border-radius: 5px;
            }
            .QRFactorizationApp QLabel {
                background: #333333;
                color: #FFFFFF;
                font-size: 13px;
                border: 1px solid #444444;
                border-radius: 4px;
                padding: 4px;
                qproperty-alignment: AlignCenter;
            }
        """)
        self.initialiser_matrice()

    def initialiser_matrice(self) -> None:
        self.layout().addWidget(self._creer_boutons_radio_decomposition())
        self.layout().addWidget(self._creer_boutons_radio_ordre())

        self._matrice_layout = QGridLayout()
        self.matrice_groupbox = self._creer_groupbox("Matrice", self._matrice_layout)
        self.layout().addWidget(self.matrice_groupbox)

        hbox_layout = QHBoxLayout()
        hbox_layout.addWidget(self._creer_boutton("Nombres entiers aléatoires", \
                                                  self.valeurs_aleatoires, 16))

        label_compris_entre = QLabel("compris entre")
        label_compris_entre.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hbox_layout.addWidget(label_compris_entre)

        self.min_spinbox = QSpinBox()
        self.min_spinbox.setRange(-99, 0)
        self.min_spinbox.setValue(-9)
        self.min_spinbox.setFixedWidth(80)
        hbox_layout.addWidget(self.min_spinbox)

        label_et = QLabel("et")
        label_et.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hbox_layout.addWidget(label_et)

        self.max_spinbox = QSpinBox()
        self.max_spinbox.setRange(0, 99)
        self.max_spinbox.setValue(9)
        self.max_spinbox.setFixedWidth(80)
        hbox_layout.addWidget(self.max_spinbox)

        hbox_layout.addWidget(self._creer_boutton("Effacer", self.effacer_entrees, 16))
        hbox_layout.addStretch()
        self.layout().addLayout(hbox_layout)

        self._q_layout = QGridLayout()
        self.matrice_q_groupbox = self._creer_groupbox("Q", self._q_layout)
        self._r_layout = QGridLayout()
        self.r_groupbox = self._creer_groupbox("R", self._r_layout)

        self._l_layout = QGridLayout()
        self.l_groupbox = self._creer_groupbox("L", self._l_layout)
        self._u_layout = QGridLayout()
        self.u_groupbox = self._creer_groupbox("U", self._u_layout)

        qr_matrices_layout = QHBoxLayout()
        qr_matrices_layout.addWidget(self.matrice_q_groupbox)
        qr_matrices_layout.addWidget(self.r_groupbox)
        self.layout().addLayout(qr_matrices_layout)

        ligne_separatrice = QFrame()
        ligne_separatrice.setFrameShape(QFrame.Shape.HLine)
        ligne_separatrice.setFrameShadow(QFrame.Shadow.Raised)
        self.layout().addWidget(ligne_separatrice)

        lu_matrices_layout = QHBoxLayout()
        lu_matrices_layout.addWidget(self.l_groupbox)
        lu_matrices_layout.addWidget(self.u_groupbox)
        self.layout().addLayout(lu_matrices_layout)

        self.matrice_q_groupbox.setVisible(self.decomposition == "QR")
        self.r_groupbox.setVisible(self.decomposition == "QR")
        self.l_groupbox.setVisible(self.decomposition == "LU")
        self.u_groupbox.setVisible(self.decomposition == "LU")

        self.layout().addWidget(self._creer_boutton("Lancer le calcul", \
                                                    self.calculer_decomposition, 16, True))
        self.layout().addWidget(self._creer_boutton("Sauvegarder le résultat", \
                                                    self.sauvegarder_resultats, 16, True))

        self.maj_entrees()
        self.maj_sorties()

    def _creer_boutons_radio_ordre(self) -> QGroupBox:
        box = QGroupBox("Sélectionner l'ordre :")
        layout = QHBoxLayout()
        self.groupe_boutons_ordre = QButtonGroup()

        for n in (2, 3, 4):
            bouton_radio = QRadioButton(f"Ordre {n}")
            bouton_radio.setProperty("ordre", n)
            bouton_radio.toggled.connect(self.changer_ordre)
            self.groupe_boutons_ordre.addButton(bouton_radio)
            layout.addWidget(bouton_radio)
            if n == self.ordre:
                bouton_radio.setChecked(True)
        layout.addStretch()
        box.setLayout(layout)
        return box

    def _creer_boutons_radio_decomposition(self) -> QGroupBox:
        box = QGroupBox("Sélectionner la décomposition :")
        layout = QHBoxLayout()
        self.groupe_boutons_decomposition = QButtonGroup()

        boutons_radio_qr = QRadioButton("QR")
        boutons_radio_qr.setProperty("type", "QR")
        boutons_radio_qr.toggled.connect(self.changer_decomposition)
        self.groupe_boutons_decomposition.addButton(boutons_radio_qr)
        layout.addWidget(boutons_radio_qr)
        if self.decomposition == "QR":
            boutons_radio_qr.setChecked(True)

        boutons_radio_lu = QRadioButton("LU (Crout)")
        boutons_radio_lu.setProperty("type", "LU")
        boutons_radio_lu.toggled.connect(self.changer_decomposition)
        self.groupe_boutons_decomposition.addButton(boutons_radio_lu)
        layout.addWidget(boutons_radio_lu)
        if self.decomposition == "LU":
            boutons_radio_lu.setChecked(True)

        layout.addStretch()
        box.setLayout(layout)
        return box

    def _creer_groupbox(self, titre: str, layout: QGridLayout) -> QGroupBox:
        box = QGroupBox(titre)
        box.setLayout(layout)
        box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        return box

    def _creer_boutton(
            self,
            libelle: str,
            callback,
            size: int = 10,
            bold: bool = False
            ) -> QPushButton:
        bouton_push = QPushButton(libelle)
        bouton_push.clicked.connect(callback)
        font = QFont()
        font.setPointSize(size)
        font.setBold(bold)
        bouton_push.setFont(font)
        return bouton_push

    def changer_ordre(self) -> None:
        bouton = self.groupe_boutons_ordre.checkedButton()
        if bouton and (ordre := bouton.property("ordre")) != self.ordre:
            self.ordre = ordre
            self.maj_entrees()
            self.maj_sorties()

    def changer_decomposition(self) -> None:
        bouton = self.groupe_boutons_decomposition.checkedButton()
        nouveau_type = bouton.property("type")
        if bouton and nouveau_type != self.decomposition:
            self.decomposition = nouveau_type
            self.matrice_q_groupbox.setVisible(self.decomposition == "QR")
            self.r_groupbox.setVisible(self.decomposition == "QR")
            self.l_groupbox.setVisible(self.decomposition == "LU")
            self.u_groupbox.setVisible(self.decomposition == "LU")
            self.maj_sorties()

    def maj_entrees(self) -> None:
        self._effacer_layout(self._matrice_layout)
        self.entrees.clear()
        validateur_regex = QRegularExpressionValidator(QRegularExpression(r"^[-]?\d*(\.\d*)?$"))
        for i in range(self.ordre):
            ligne = []
            for j in range(self.ordre):
                case_lineedit = QLineEdit("0")
                case_lineedit.setFixedSize(75, 25)
                case_lineedit.setAlignment(Qt.AlignmentFlag.AlignCenter)
                case_lineedit.setValidator(validateur_regex)
                self._matrice_layout.addWidget(case_lineedit, i, j)
                ligne.append(case_lineedit)
            self.entrees.append(ligne)
        self.matrice_groupbox.setTitle(f"Matrice d'entrée ({self.ordre}x{self.ordre})")

    def maj_sorties(self) -> None:
        output_data = (
            (self._q_layout, self.libelles_q, "Q (orthogonale)"),
            (self._r_layout, self.libelles_r, "R (triangulaire supérieure)"),
            (self._l_layout, self.libelles_l, "L (triangulaire inférieure unitaire)"),
            (self._u_layout, self.libelles_u, "U (triangulaire supérieure)")
        )

        for layout, libelles, title_suffix in output_data:
            self._effacer_layout(layout)
            libelles.clear()
            for i in range(self.ordre):
                ligne = []
                for j in range(self.ordre):
                    libelle = QLabel("-")
                    libelle.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    libelle.setMinimumSize(72, 28)
                    libelle.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
                    layout.addWidget(libelle, i, j)
                    ligne.append(libelle)
                libelles.append(ligne)
            if layout == self._q_layout:
                self.matrice_q_groupbox.setTitle(f"{title_suffix} ({self.ordre}x{self.ordre})")
            elif layout == self._r_layout:
                self.r_groupbox.setTitle(f"{title_suffix} ({self.ordre}x{self.ordre})")
            elif layout == self._l_layout:
                self.l_groupbox.setTitle(f"{title_suffix} ({self.ordre}x{self.ordre})")
            elif layout == self._u_layout:
                self.u_groupbox.setTitle(f"{title_suffix} ({self.ordre}x{self.ordre})")

    def _effacer_layout(self, layout) -> None:
        while layout.count():
            enfant = layout.takeAt(0)
            if enfant.widget():
                enfant.widget().deleteLater()

    def calculer_decomposition(self) -> None:
        try:
            matrice = np.array([[float(case.text()) for case in ligne] for ligne in self.entrees])
        except ValueError:
            QMessageBox.warning(self, "Erreur", "Entrée invalide dans la matrice.")
            return

        if self.decomposition == "QR" and np.isclose(np.linalg.det(matrice), 0):
            QMessageBox.warning(self, "Erreur", "La matrice est singulière pour la méthode QR.")
            return

        if self.decomposition == "LU" and np.isclose(np.linalg.det(matrice), 0):
            QMessageBox.warning(self, "Erreur", "La matrice est singulière pour la méthode LU.")
            return

        for i in range(self.ordre):
            for j in range(self.ordre):
                self.libelles_q[i][j].setText("-")
                self.libelles_r[i][j].setText("-")
                self.libelles_l[i][j].setText("-")
                self.libelles_u[i][j].setText("-")

        if self.decomposition == "QR":
            try:
                Q, R = spla.qr(matrice)
                for i in range(self.ordre):
                    for j in range(self.ordre):
                        self.libelles_q[i][j].setText(f"{Q[i, j]:.4f}")
                        self.libelles_r[i][j].setText(f"{R[i, j]:.4f}")
            except spla.LinAlgError as e:
                QMessageBox.warning(self, "Erreur QR", f"Décomposition QR impossible : {e}")
                return
        elif self.decomposition == "LU":
            try:
                _, L, U = spla.lu(matrice) # on ne s'occupe pas de la matrice de permutation ici
                for i in range(self.ordre):
                    for j in range(self.ordre):
                        self.libelles_l[i][j].setText(f"{L[i, j]:.4f}")
                        self.libelles_u[i][j].setText(f"{U[i, j]:.4f}")
            except spla.LinAlgError as e:
                QMessageBox.warning(self, "Erreur LU", f"Décomposition LU impossible : {e}")
                return

    def sauvegarder_resultats(self) -> None:
        matrices_a_sauvegarder = {}
        if self.decomposition == "QR":
            matrices_a_sauvegarder["Q"] = self._recuperer_matrice_depuis_libelles(self.libelles_q)
            matrices_a_sauvegarder["R"] = self._recuperer_matrice_depuis_libelles(self.libelles_r)
        elif self.decomposition == "LU":
            matrices_a_sauvegarder["L"] = self._recuperer_matrice_depuis_libelles(self.libelles_l)
            matrices_a_sauvegarder["U"] = self._recuperer_matrice_depuis_libelles(self.libelles_u)

        if all(matrix is None for matrix in matrices_a_sauvegarder.values()):
            QMessageBox.information(self, "Aucune donnée", "Veuillez d'abord lancer un calcul \
                                    avant de sauvegarder.")
            return

        nom_fichier, _ = QFileDialog.getSaveFileName(self, "Sauvegarder le résultat (.csv)", \
                                                     "résultats.csv", "Fichiers CSV (*.csv);;Tous fichiers (*)")
        if nom_fichier:
            try:
                with open(nom_fichier, "w", newline="", encoding="utf-8") as fichier:
                    for nom, matrice in matrices_a_sauvegarder.items():
                        if matrice is not None:
                            fichier.write(f"--> Matrice {nom}\n")
                            np.savetxt(fichier, matrice, fmt="%.4f", delimiter=";", \
                                       newline="\n", comments="")
                            fichier.write("\n")
                QMessageBox.information(self, "Sauvegarde réussie", f"Les résultats ont été \
                                        sauvegardés dans :\n{nom_fichier}")
            except Exception as e:
                QMessageBox.critical(self, "Erreur de sauvegarde", f"Une erreur est survenue \
                                     lors de la sauvegarde du fichier :\n{e}")

    def _recuperer_matrice_depuis_libelles(self, liste_libelles: list) -> np.ndarray | None:
        if not liste_libelles or not liste_libelles[0]:
            return None
        
        if liste_libelles[0][0].text() == "-":
            return None

        matrix = np.zeros((self.ordre, self.ordre))
        try:
            for i in range(self.ordre):
                for j in range(self.ordre):
                    matrix[i, j] = float(liste_libelles[i][j].text())
            return matrix
        except ValueError:
            return None

    def valeurs_aleatoires(self) -> None:
        rng = np.random.default_rng()
        for ligne in self.entrees:
            for case in ligne:
                case.setText(str(rng.integers(self.min_spinbox.value(), self.max_spinbox.value())))

    def effacer_entrees(self) -> None:
        for ligne in self.entrees:
            for case in ligne:
                case.setText("0")
        for i in range(self.ordre):
            for j in range(self.ordre):
                self.libelles_q[i][j].setText("-")
                self.libelles_r[i][j].setText("-")
                self.libelles_l[i][j].setText("-")
                self.libelles_u[i][j].setText("-")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = FactorisationQRLU()
    main.show()
    sys.exit(app.exec())