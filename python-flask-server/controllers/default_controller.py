from rdkit import Chem
import rdkit.Chem.inchi
from rdkit.Chem.Descriptors import MolWt
from urllib.error import HTTPError
import cirpy
import re
from werkzeug.contrib.cache import SimpleCache
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

cas_to_inchi_cache = SimpleCache()
cas_to_inchikey_cache = SimpleCache()
cas_to_smiles_cache = SimpleCache()
cas_to_names_cache = SimpleCache()
inchi_to_cas_cache = SimpleCache()
inchi_to_names_cache = SimpleCache()
inchikey_to_cas_cache = SimpleCache()
inchikey_to_smiles_cache = SimpleCache()
inchikey_to_inchi_cache = SimpleCache()
inchikey_to_names_cache = SimpleCache()
smiles_to_cas_cache = SimpleCache()
smiles_to_names_cache = SimpleCache()
name_to_inchi_cache = SimpleCache()
name_to_inchikey_cache = SimpleCache()
name_to_smiles_cache = SimpleCache()
name_to_cas_cache = SimpleCache()


caches = {}
caches['cas'] = {}
caches['cas']['smiles'] = cas_to_smiles_cache
caches['cas']['stdinchi'] = cas_to_inchi_cache
caches['cas']['stdinchikey'] = cas_to_inchikey_cache
caches['cas']['names'] = cas_to_names_cache
caches['stdinchi'] = {}
caches['stdinchi']['cas'] = inchi_to_cas_cache
caches['stdinchi']['names'] = inchi_to_names_cache
caches['stdinchikey'] = {}
caches['stdinchikey']['cas'] = inchikey_to_cas_cache
caches['stdinchikey']['stdinchi'] = inchikey_to_inchi_cache
caches['stdinchikey']['smiles'] = inchikey_to_smiles_cache
caches['stdinchikey']['names'] = inchikey_to_names_cache
caches['smiles'] = {}
caches['smiles']['cas'] = smiles_to_cas_cache
caches['smiles']['names'] = smiles_to_names_cache
caches['name'] = {}
caches['name']['smiles'] = name_to_smiles_cache
caches['name']['stdinchi'] = name_to_inchi_cache
caches['name']['stdinchikey'] = name_to_inchikey_cache
caches['name']['cas'] = name_to_cas_cache


# In-Process conversions using rdkit:
def as_svg_get(smiles, width, height):
    mol = Chem.MolFromSmiles(smiles)
    if (mol is None):
        return ("Could not parse input: " + smiles, 500)

    kekulize = True
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return svg.replace('svg:', '').replace('xmlns:svg=', 'xmlns=')


def mol_weight_get(smiles):
    m = Chem.MolFromSmiles(smiles)
    if (m is None):
        return ("Could not parse input: " + smiles, 500)

    return {"smiles": smiles, "molWeight": MolWt(m)}


def inchi_to_inchikey_get(inchi):
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    if inchikey is None:
        return ("Could not parse input: " + inchi, 500)

    return {"inchikey": inchikey}


def inchi_to_smiles_get(inchi):
    m = rdkit.Chem.inchi.MolFromInchi(inchi)
    if (m is None):
        return ("Could not parse input: " + inchi, 500)

    return {"smiles": Chem.MolToSmiles(m)}


def smiles_to_inchi_get(smiles):
    m = Chem.MolFromSmiles(smiles)
    if (m is None):
        return ("Could not parse input: " + smiles, 500)

    return {"inchi": rdkit.Chem.inchi.MolToInchi(m)}


def smiles_to_inchikey_get(smiles):
    m = Chem.MolFromSmiles(smiles)
    if (m is None):
        return ("Could not parse input: " + smiles, 500)

    inchi = rdkit.Chem.inchi.MolToInchi(m)
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    return {"inchikey": inchikey}

# conversions using external REST services


class CirpyError(Exception):
    def __init__(self, code, message):
        self.code = code
        self.message = message


def resolve_via_cirpy(identifier, target, source):
    try:
        converted = caches[source][target].get(identifier)
        if converted is None:
            sourcehint = 'cas_number' if source == 'cas' else source
            sourcehints = ['name_by_opsin', 'name_by_cir'] if sourcehint == 'name' else [sourcehint]
            converted = cirpy.resolve(identifier, target, sourcehints)
            caches[source][target].set(identifier, converted)
        return converted
    except HTTPError as err:
        if err.code == 504 or err.code == 408:
            raise CirpyError(504, "Timeout while waiting for identifier resolution service")
        raise CirpyError(500, "HTTPError while communicating with identifier resolution service" + err.reason)


def inchikey_to_inchi_get(inchikey):
    try:
        inchi = resolve_via_cirpy(inchikey, 'stdinchi', 'stdinchikey')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_names_get(inchikey):
    try:
        names = resolve_via_cirpy(inchikey, 'names', 'stdinchikey')
        return {"names": names}
    except CirpyError as err:
        return (err.message, err.code)



def inchi_to_cas_get(inchi):
    try:
        cas = resolve_via_cirpy(inchi, 'cas', 'stdinchi')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def inchi_to_names_get(inchi):
    try:
        names = resolve_via_cirpy(inchi, 'names', 'stdinchi')
        return {"names": names}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_cas_get(inchikey):
    try:
        cas = resolve_via_cirpy(inchikey, 'cas', 'stdinchikey')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchi_get(cas):
    try:
        inchi = resolve_via_cirpy(cas, 'stdinchi', 'cas')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchikey_get(cas):
    try:
        inchikey = clean_inchi_key(resolve_via_cirpy(cas, 'stdinchikey', 'cas'))
        return {"inchikey": inchikey}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_smiles_get(cas):
    try:
        smiles = resolve_via_cirpy(cas, 'smiles', 'cas')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_names_get(cas):
    try:
        names = resolve_via_cirpy(cas, 'names', 'cas')
        return {"names": names}
    except CirpyError as err:
        return (err.message, err.code)


def smiles_to_cas_get(smiles):
    try:
        cas = resolve_via_cirpy(smiles, 'cas', 'smiles')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def smiles_to_names_get(smiles):
    try:
        names = resolve_via_cirpy(smiles, 'names', 'smiles')
        return {"names": names}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_smiles_get(inchikey):
    try:
        smiles = resolve_via_cirpy(inchikey, 'smiles', 'stdinchikey')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)


def name_to_smiles_get(name):
    try:
        smiles = resolve_via_cirpy(name, 'smiles', 'name')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)


def name_to_inchi_get(name):
    try:
        inchi = resolve_via_cirpy(name, 'stdinchi', 'name')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)


def name_to_inchikey_get(name):
    try:
        inchikey = resolve_via_cirpy(name, 'stdinchikey', 'name')
        return {"inchikey": inchikey}
    except CirpyError as err:
        return (err.message, err.code)


def name_to_cas_get(name):
    try:
        cas = resolve_via_cirpy(name, 'cas', 'name')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def clean_inchi_key(inchikey):
    return re.sub("(?i)InChiKey=", "", inchikey)
