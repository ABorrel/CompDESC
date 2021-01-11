from rdkit.Chem import rdMolDescriptors




def getMQN(mol):
    l_MQNs = rdMolDescriptors.MQNs_(mol)
    return l_MQNs


_mqn = {"MQN1": getMQN,
        "MQN2": getMQN,
        "MQN3": getMQN,
        "MQN4": getMQN,
        "MQN5": getMQN,
        "MQN6": getMQN,
        "MQN7": getMQN,
        "MQN8": getMQN,
        "MQN9": getMQN,
        "MQN10": getMQN,
        "MQN11": getMQN,
        "MQN12": getMQN,
        "MQN13": getMQN,
        "MQN14": getMQN,
        "MQN15": getMQN,
        "MQN16": getMQN,
        "MQN17": getMQN,
        "MQN18": getMQN,
        "MQN19": getMQN,
        "MQN20": getMQN,
        "MQN21": getMQN,
        "MQN22": getMQN,
        "MQN23": getMQN,
        "MQN24": getMQN,
        "MQN25": getMQN,
        "MQN26": getMQN,
        "MQN27": getMQN,
        "MQN28": getMQN,
        "MQN29": getMQN,
        "MQN30": getMQN,
        "MQN31": getMQN,
        "MQN32": getMQN,
        "MQN33": getMQN,
        "MQN34": getMQN,
        "MQN35": getMQN,
        "MQN36": getMQN,
        "MQN37": getMQN,
        "MQN38": getMQN,
        "MQN39": getMQN,
        "MQN40": getMQN,
        "MQN41": getMQN,
        "MQN42": getMQN}


def GetMQNs(mol):
    
    l_MQNs = getMQN(mol)
    d_out = {}
    i = 0
    imax = len(l_MQNs)
    while i < imax:
        d_out["MQN%s"%(i+1)] = l_MQNs[i]
        i = i + 1

    return d_out
    