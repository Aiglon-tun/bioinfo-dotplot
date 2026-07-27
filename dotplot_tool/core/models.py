from dataclasses import dataclass 

@dataclass(frozen=True) #Ce decorateur rend la classe immuable, cad que les instances de cette classe ne peuvent
#pas être modifiées après leur création.


class GenomeCandidate :
    assembly_accession : str
    organism_name : str
    assembly_level : str = ""
    refseq_category : str = ""