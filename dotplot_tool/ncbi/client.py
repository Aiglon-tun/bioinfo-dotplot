from __future__ import annotations

import time
import requests


class NcbiClient:
    """
    Client minimal pour interroger les eUtils du NCBI.
    """

    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(
        self,
        email: str = "student@example.com",
        api_key: str | None = None,
        timeout: int = 30,
    ) -> None:
        self.email = email
        self.api_key = api_key
        self.timeout = timeout
        self.session = requests.Session()

    def _get(self, path: str, params: dict) -> dict:
        params = dict(params)
        params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{self.BASE}/{path}"

        for attempt in range(3):
            r = self.session.get(url, params=params, timeout=self.timeout)
            if r.status_code == 200:
                return r.json()
            time.sleep(0.7 * (attempt + 1))

        r.raise_for_status()
        raise RuntimeError("Erreur lors de la requête NCBI")

    def esearch(self, db: str, term: str, retmax: int = 20) -> list[str]:
        data = self._get(
            "esearch.fcgi",
            {
                "db": db,
                "term": term,
                "retmode": "json",
                "retmax": str(retmax),
                "sort": "relevance",
            },
        )
        return data.get("esearchresult", {}).get("idlist", [])

    def esummary(self, db: str, ids: list[str]) -> dict:
        if not ids:
            return {}

        data = self._get(
            "esummary.fcgi",
            {
                "db": db,
                "id": ",".join(ids),
                "retmode": "json",
            },
        )
        return data.get("result", {})
