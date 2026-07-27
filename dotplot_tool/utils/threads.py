from __future__ import annotations
from PySide6.QtCore import QObject, Signal, QRunnable, Slot


class WorkerSignals(QObject):
    """Signaux Qt émis par un worker exécuté en arrière-plan.

    On sépare les signaux dans un QObject dédié car QRunnable n'est pas
    un QObject. Cela permet de connecter facilement ces signaux à l'UI.
    """

    # Émis quand la tâche est totalement terminée (succès ou erreur).
    finished = Signal()
    # Émis en cas d'exception: on envoie un message texte simple.
    error = Signal(str)
    # Émis quand la fonction a produit un résultat.
    result = Signal(object)


class Worker(QRunnable):
    """Wrapper générique pour exécuter une fonction dans le QThreadPool.

    L'idée: prendre n'importe quelle fonction Python + ses arguments,
    l'exécuter hors du thread UI, puis renvoyer résultat/erreur via signaux.
    """

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        # Fonction métier à exécuter (ex: requête réseau NCBI).
        self.fn = fn
        # Arguments positionnels transmis à la fonction.
        self.args = args
        # Arguments nommés transmis à la fonction.
        self.kwargs = kwargs
        # Objet de signaux auquel l'interface peut se connecter.
        self.signals = WorkerSignals()

    @Slot()
    def run(self):
        # Méthode appelée par Qt quand le worker démarre dans le thread pool.
        try:
            out = self.fn(*self.args, **self.kwargs)
            self.signals.result.emit(out)
        except Exception as e:
            # On remonte l'erreur à l'UI sans faire crasher l'application.
            self.signals.error.emit(str(e))
        finally:
            # Toujours émis pour permettre le nettoyage UI (boutons, loaders...).
            self.signals.finished.emit()
