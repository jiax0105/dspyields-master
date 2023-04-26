import sys
sys.path.append('.')

import json
import os
import csv
import pandas as pd
import numpy as np
import yaml


def write_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)


def read_json(file):
    with open(file, 'r') as f:
        data = json.load(f)
    return data


# def get_config(yaml_file='./config/config.yaml'):
#     with open(yaml_file, 'r', encoding='utf-8') as f:
#         data_dict = yaml.load(f, Loader=yaml.FullLoader)
#     return data_dict
def get_config():
    #os.path.dirname(__file__) #当前文件的路径
    #os.path.dirname(os.path.dirname(__file__)) #当前文件的上一级目录
    project_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))  # 上级目录
    yaml_file = os.path.join(project_dir, 'config/config.yaml')
    with open(yaml_file, 'r', encoding='utf-8') as f:
        data_dict = yaml.load(f, Loader=yaml.FullLoader)

    for key in data_dict['file_path']:
        data_dict['file_path'][key] = os.path.join(project_dir, data_dict['file_path'][key])

    return data_dict

class SingleReaction(object):
    def __init__(self, S, P, reaction, Tanimoto=0.0, next=None):
        self.S = S
        self.P = P
        self.reaction = reaction
        self.Tanimoto = Tanimoto
        self.next = next

class SingleLinkList(object):
    def __init__(self):
        self._head = None

    def is_empty(self):
        return self._head == None

    def clear(self):
        self._head = None

    def get_rxn_dict(self):
        rxn_dict = {}
        cur = self._head
        while cur != None:
            rxn_dict[cur.reaction['kegg_id']] = cur.reaction['equation']
            cur = cur.next
        return rxn_dict

    def get_rxn_list(self):
        rxn_list = []
        cur = self._head
        while cur != None:
            rxn_list.append(cur.reaction['kegg_id'])
            cur = cur.next
        return rxn_list

    def get_cpd_list(self):
        cpd_list = []
        cur = self._head
        while cur != None:
            if cur.S not in cpd_list:
                cpd_list.append(cur.S)
            if cur.P not in cpd_list:
                cpd_list.append(cur.P)
            cur = cur.next
        return cpd_list

    def travel(self):
        """遍历整个链表"""
        travel_list = []
        cur = self._head
        while cur != None:
            travel_list.append([cur.reaction['kegg_id'], cur.S, cur.P])
            print(cur.reaction['kegg_id'], ":", cur.S, " --> ", cur.P, ",")
            cur = cur.next
        return travel_list


    def length(self):
        count = 0
        cur = self._head
        while cur != None:
            count += 1
            cur = cur.next
        return count

    def add(self, node):
        """在头部添加节点，节点是已知，更新头节点"""
        # node = SingleNode()
        # node = Node(item)
        node.next = self._head
        self._head = node

    def get_first_node(self):
        return self._head

    def append(self, node):
        """在尾部添加节点"""
        if self.is_empty():
            self._head = node
        else:
            cur = self._head
            # node = Node(item)
            while cur.next != None:
                cur = cur.next
            cur.next = node

    def insert(self, pos, node):
        """在选定的位置添加节点"""
        cur = self._head
        # node = Node(item)
        count = 0
        if pos <= 0:
            self.add(node)
        elif pos > (self.length() - 1):
            self.append(node)
        else:
            while count < (pos - 1):
                count += 1
                cur = cur.next
            node.next = cur.next
            cur.next = node


def main():
    data = get_config()
    print(data)

if __name__ == "__main__":
    main()