import re
import sys
import time
import queue
import threading
from time import sleep
from bs4 import BeautifulSoup
import requests
from requests import get as rget
from rdkit import Chem
import requests
import numpy as np


# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def kegg_get(kegg_id, server="http://rest.kegg.jp"):
    """
    Downloads the raw REST text entry for the provided KEGG ID,
    via the KEGG rest API @ https://rest.kegg.jp/
    """
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36'
    }
    # Construct query
    rest_address = "/".join([server, "get", kegg_id])
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        session = requests.Session()
        r = session.get(rest_address, headers=headers)
        # r = rget(rest_address, headers=headers)
        if r.status_code == 200:
            return r.text
        else:
            # Server not returning a result, try again
            time.sleep(2)
    # The connection attempt limit was reached
    sError("Warning: Unable to download KEGG data for '%s'.\n" % str(kegg_id))
    return None


def threaded_kegg_get(queries):
    """Threaded implementation of kegg_get."""

    def worker():
        while True:
            query = work.get()
            if query is None:
                break
            result = kegg_get(query)
            output.put((query, result))
            work.task_done()

    # Initialise queues
    work = queue.Queue()
    output = queue.Queue()

    for query in queries:
        work.put(query)

    # Start threads
    threads = []
    for i in range(16):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # Report on progress
    while True:
        if len(queries) == 0:
            progress = 100.0
        else:
            progress = float(output.qsize() / len(queries) * 100)
        sys.stdout.write("\rQuerying KEGG... %0.1f%%" % progress)
        sys.stdout.flush()
        if output.qsize() == len(queries):
            print("")
            break
        sleep(0.5)

    # Join work queue
    work.join()

    # Stop workers
    for i in range(len(threads)):
        work.put(None)
    for t in threads:
        t.join()

    # Get results
    results = {}
    while not output.empty():
        result = output.get()
        results[result[0]] = result[1]
    return results


def create_kegg_dict(kegg_text):
    """
    Parses a KEGG REST text record into a dictionary. Accepts a single record,
    as it stops after the first '///'. All lines under a key are presented as a
    list of lists.
    """
    # Initialize dictionary
    kegg_dict = {}
    # Iterate over lines in raw text
    for line in kegg_text.split("\n"):
        if line == "///":
            # The first entry ends here
            break
        line = re.split(" +", line.rstrip())
        if line[0] != "":
            key = line[0]
            line = line[1:]
            try:
                kegg_dict[key].append(line)
            except KeyError:
                kegg_dict[key] = [line]
        else:
            kegg_dict[key].append(line[1:])
    return kegg_dict


def kegg_smiles(kegg_id, server="http://rest.kegg.jp"):
    """Downloads a KEGG compound molecule object and converts it to SMILES."""
    # Ensure that the ID is a valid KEGG compound
    if not re.fullmatch("^C[0-9]{5}$", kegg_id):
        sError("\nWarning: '%s' is not a valid KEGG compound ID.\n" % str(kegg_id))
        return None
    # Set up the query
    rest_address = "/".join([server, "get", "cpd:" + kegg_id, "mol"])
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36'}
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        r = rget(rest_address, headers=headers)
        if r.status_code == 200:
            try:
                mol = Chem.MolFromMolBlock(r.text)
                smiles = Chem.MolToSmiles(mol)
                return smiles
            except:
                sError("\nWarning: SMILES could not be produced for KEGG ID '%s'.\n" % str(kegg_id))
                return None
        else:
            time.sleep(2)
    # The connection attempt limit was reached
    sError("\nWarning: Unable to download SMILES of molecule data for '%s'.\n" % str(kegg_id))
    return None


def KeggId2BiggId(kegg_id):
    res = re.search('R|C', kegg_id).group()

    if res == 'R':
        cookie = '_ga=GA1.2.208736496.1642561994; _gid=GA1.2.426276517.1656811686; _gat=1'
        database_source = 'kegg.reaction'
        re_pat = r'<a href="/models/universal/reactions/.*">(.*)</a>'

    elif res == 'C':
        cookie = '_ga=GA1.2.208736496.1642561994; _gid=GA1.2.842153841.1656397739; _gat=1'
        database_source = 'kegg.compound'
        re_pat = r'<a href="/models/universal/metabolites/.*">(.*)</a>'
    else:
        print('%s is a fault KEGG ID' % (kegg_id))
        return

    url = 'http://bigg.ucsd.edu/advanced_search_external_id_results'

    headers = {
        'Cookie': cookie,
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36'
    }

    payload = {
        'database_source': database_source,
        'query': kegg_id
    }

    con_attempts = 0
    while con_attempts < 5:
        con_attempts += 1
        r = requests.post(url=url, data=payload, headers=headers)
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, features="html.parser")
            lines = soup.find_all('a', href=True)[5:-1]
            result = [re.findall(re_pat, str(line)) for line in lines]
            if result == []:
                return None
            else:
                bigg_ids = [res[0] for res in result]
                return bigg_ids

        else:
            time.sleep(2)

    return None


class MyReaction():
    """
    Basic KEGG reaction object constructed from a reaction dictionary.

    equation        The KEGG equation string
    compound_ids    A set of KEGG IDs representing reactants and products
    reactants       A set of KEGG IDs representing reactants only
    products        A set of KEGG IDs representing products only
    """

    def __init__(self, kegg_id):
        self.kegg_dict = create_kegg_dict(kegg_get(kegg_id))

        # Reaction ID
        self.kegg_id = kegg_id  # self.kegg_dict["ENTRY"][0][0]

        self.bigg_id = KeggId2BiggId(self.kegg_id)

        # Get the equation list
        eq_list = self.kegg_dict["EQUATION"][0]
        # Original equation string
        self.equation = " ".join(eq_list)

        # Get the definition list
        def_list = self.kegg_dict["DEFINITION"][0]
        self.definition = " ".join(def_list)

        # Create sets of compounds, reactants and products
        is_comp = re.compile("^C[0-9]{5}")
        self.compounds = set(filter(is_comp.match, eq_list))
        self.reactants = set(filter(is_comp.match, eq_list[0:eq_list.index("<=>")]))
        self.products = set(filter(is_comp.match, eq_list[eq_list.index("<=>") + 1:]))

        # Get Enzyme of reaction
        # self.enzyme = self.kegg_dict["ENZYME"][0]


class MyCompound():

    def __init__(self, kegg_id) -> None:

        self.kegg_dict = create_kegg_dict(kegg_get(kegg_id))

        self.kegg_id = kegg_id

        self.bigg_id = KeggId2BiggId(self.kegg_id)

        self.formula = self._get_formula()

        self.name = self._get_name()

        self.cid = self.getPubchemCid()

        self.chebi = self.getChebi()

        # self.smile = self.get_smile()

    def _get_formula(self):

        if "FORMULA" in self.kegg_dict.keys():
            return self.kegg_dict["FORMULA"][0][0]

        else:
            return None

    def _get_name(self):

        names = []

        for i in range(len(self.kegg_dict["NAME"])):

            if len(self.kegg_dict["NAME"][i]) > 1:

                name = " ".join(self.kegg_dict["NAME"][i]).strip(";")

            else:

                name = self.kegg_dict["NAME"][i][0].strip(";")

            names.append(name)

        return names

    def getPubchemCid(self):

        if "DBLINKS" in self.kegg_dict.keys():
            for ids in self.kegg_dict["DBLINKS"]:
                if ids[0] == 'PubChem:': return int(ids[1])
        return None

    def getChebi(self):

        if "DBLINKS" in self.kegg_dict.keys():
            for ids in self.kegg_dict["DBLINKS"]:
                if ids[0] == 'ChEBI:': return int(ids[1])
        return None

    def get_smile(self):
        return kegg_smiles(self.kegg_id)


def parse_equation(equation):
    """Parses a KEGG reaction string into a tuple of lists of tuples."""
    # Split the equation into lists
    eq_list = re.split(" +", equation)
    reactants = eq_list[0:eq_list.index("<=>")]
    products = eq_list[eq_list.index("<=>") + 1:]

    # Set up parser
    def stoichiometry_parse(stoichiometry_list):
        output = []
        s = 1
        for segment in stoichiometry_list:
            if re.match("^[0-9]+$", segment):
                s = int(segment)
                continue
            if segment == "+":
                continue
            else:
                output.append((s, segment))
                s = 1
        return output

    return (stoichiometry_parse(reactants), stoichiometry_parse(products))


def parser_equation_coef(equation):
    eq_dict = {}
    S, P = equation.rstrip().split(' <=> ')
    S = S.split(' + ')
    pat = '^\d+\ |^[1-9nm]+\ |^\([1-9nm]+\)\ |\([1-9nmx]+\)$|^\([nm][+-][nm\d]\)\ |^[nm][+-][nm1-3]\ |\([nm][+-][nmx\d]\)$'
    for s in S[:]:
        coef = re.findall(pat, s)
        # REACOMODO COEFICIENTES
        if len(coef) > 0:
            s = s.strip().replace(coef[0], '', 1)
        else:
            coef = ['1']
        # ----------------------------------------------------
        # CORRECCIONES 更正
        # ----------------------
        coef = coef[0].strip()
        coef = coef.replace('(n-2)', '98')
        coef = coef.replace('(n-1)', '99')
        coef = coef.replace('(n)', '100')
        coef = coef.replace('(n+1)', '101')
        coef = coef.replace('(n+2)', '102')
        coef = coef.replace('(m-2)', '98')
        coef = coef.replace('(m-1)', '99')
        coef = coef.replace('(m)', '100')
        coef = coef.replace('(m+1)', '101')
        coef = coef.replace('(m+2)', '102')
        coef = coef.replace('(n+m)', '200')
        coef = coef.replace('(m+n)', '200')
        coef = coef.replace('(n-x)', '90')
        coef = coef.replace('n-1', '99')
        coef = coef.replace('n+1', '101')
        coef = coef.replace('n', '100')
        coef = coef.replace('(x)', '10')

        eq_dict[s] = float(coef[0]) * -1

    # PROCESAR Y PARSEAR PRODUCTS
    P = P.split(' + ')
    for p in P[:]:
        coef = re.findall(pat, p)

        # REACOMODO COEFICIENTES
        if len(coef) > 0:
            p = p.replace(coef[0], '', 1)
        else:
            coef = ['1']

        coef = coef[0].strip()
        coef = coef.replace('(n-2)', '98')
        coef = coef.replace('(n-1)', '99')
        coef = coef.replace('(n)', '100')
        coef = coef.replace('(n+1)', '101')
        coef = coef.replace('(n+2)', '102')
        coef = coef.replace('(m-2)', '98')
        coef = coef.replace('(m-1)', '99')
        coef = coef.replace('(m)', '100')
        coef = coef.replace('(m+1)', '101')
        coef = coef.replace('(m+2)', '102')
        coef = coef.replace('(n+m)', '200')
        coef = coef.replace('(m+n)', '200')
        coef = coef.replace('(n-x)', '90')
        coef = coef.replace('n-1', '99')
        coef = coef.replace('n+1', '101')
        coef = coef.replace('n', '100')
        coef = coef.replace('(x)', '10')

        eq_dict[p] = float(coef[0]) * 1

    return eq_dict


# def split_equation(equation):
#     eq_list = re.split(" +", equation)
#     is_comp = re.compile("^C[0-9]{5}")
#     def check_func(string):
#         return is_comp.match(string).group()
#     reactants = list(set(filter(check_func, eq_list[0:eq_list.index("<=>")])))
#     products = list(set(filter(check_func, eq_list[eq_list.index("<=>")+1:])))
#     return reactants, products


def split_equation(equation):
    _reactants, _products = equation.split(' <=> ')
    pat = re.compile("C[0-9]{5}")
    reactants = list(set(pat.findall(_reactants)))
    products = list(set(pat.findall(_products)))
    return reactants, products


if __name__ == '__main__':
    # print(KeggId2BiggId('C00001'))

    # rxn = 'R04936'
    # r = MyReaction(rxn)

    # eq = "C00677 + C00039(n) <=> C00013 + C00039(n+1)" #R00379
    # eq = 'C00404 + n C00001 <=> (n+1) C02174' #R00001
    eq = 'C00043 + 2 C00003 + C00001 <=> C04573 + 2 C00004 + 2 C00080'  # R00421
    print(split_equation(eq))
    print(parse_equation(eq))
