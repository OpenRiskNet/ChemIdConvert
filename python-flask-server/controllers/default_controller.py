from rdkit import Chem
import rdkit.Chem.inchi
from urllib.error import HTTPError
import cirpy
import re
from werkzeug.contrib.cache import SimpleCache

cas_to_inchi_cache = SimpleCache()
cas_to_inchikey_cache = SimpleCache()
cas_to_smiles_cache = SimpleCache()
inchi_to_cas_cache = SimpleCache()
inchikey_to_cas_cache = SimpleCache()
inchikey_to_smiles_cache = SimpleCache()
inchikey_to_inchi_cache = SimpleCache()
smiles_to_cas_cache = SimpleCache()

caches = {}
caches['cas'] = {}
caches['cas']['smiles'] = cas_to_smiles_cache
caches['cas']['stdinchi'] = cas_to_inchi_cache
caches['cas']['stdinchikey'] = cas_to_inchikey_cache
caches['stdinchi'] = {}
caches['stdinchi']['cas'] = inchi_to_cas_cache
caches['stdinchikey'] = {}
caches['stdinchikey']['cas'] = inchikey_to_cas_cache
caches['stdinchikey']['stdinchi'] = inchikey_to_inchi_cache
caches['stdinchikey']['smiles'] = inchikey_to_smiles_cache
caches['smiles'] = {}
caches['smiles']['cas'] = smiles_to_cas_cache


# In-Process conversions using rdkit:

def inchi_to_inchikey_get(inchi: str) -> str:
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    if inchikey == None:
        return ("Could not parse input: " + inchi, 500)

    return {"inchikey": inchikey}


def inchi_to_smiles_get(inchi: str) -> str:
    m = rdkit.Chem.inchi.MolFromInchi(inchi)
    if (m == None):
        return ("Could not parse input: " + inchi, 500)

    return {"smiles": Chem.MolToSmiles(m)}


def smiles_to_inchi_get(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    return {"inchi": rdkit.Chem.inchi.MolToInchi(m)}


def smiles_to_inchikey_get(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    inchi = rdkit.Chem.inchi.MolToInchi(m)
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    return {"inchikey": inchikey}

# conversions using external REST services


class CirpyError(Exception):
    def __init__(self, code, message):
        self.code = code
        self.message = message


def resolve_via_cirpy(identifier: str, target: str, source: str,) -> str:
    try:
        converted = caches[source][target].get(identifier)
        if converted is None:
            sourcehint = 'cas_number' if source == 'cas' else source
            converted = cirpy.resolve(identifier, target, [sourcehint])
            caches[source][target].set(identifier, converted)
        return converted
    except HTTPError as err:
        if err.code == 504 or err.code == 408:
            raise CirpyError(504, "Timeout while waiting for identifier resolution service")
        raise CirpyError(500, "HTTPError while communicating with identifier resolution service" + err.reason)


def inchikey_to_inchi_get(inchikey) -> str:
    try:
        inchi = resolve_via_cirpy(inchikey, 'stdinchi', 'stdinchikey')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)



def inchi_to_cas_get(inchi) -> str:
    try:
        cas = resolve_via_cirpy(inchi, 'cas', 'stdinchi')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_cas_get(inchikey) -> str:
    try:
        cas = resolve_via_cirpy(inchikey, 'cas', 'stdinchikey')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchi_get(cas) -> str:
    try:
        inchi = resolve_via_cirpy(cas, 'stdinchi', 'cas')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchikey_get(cas) -> str:
    try:
        inchikey = clean_inchi_key(resolve_via_cirpy(cas, 'stdinchikey', 'cas'))
        return {"inchikey": inchikey}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_smiles_get(cas) -> str:
    try:
        smiles = resolve_via_cirpy(cas, 'smiles', 'cas')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)


def smiles_to_cas_get(smiles) -> str:
    try:
        cas = resolve_via_cirpy(smiles, 'cas', 'smiles')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_smiles_get(inchikey) -> str:
    try:
        smiles = resolve_via_cirpy(inchikey, 'smiles', 'stdinchikey')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)

def clean_inchi_key(inchikey):
    return re.sub("(?i)InChiKey=", "", inchikey)
