# coding:utf-8
import os.path
from os.path import split

import numpy as np
import networkx as nx
import sys
import matplotlib.pyplot as plt
import json

from numpy.matrixlib.defmatrix import matrix



def labeling(input_tsv, plp, blb, feature, transList):
	sites_p = {}
	sites_b = {}
	linenum = 0
	# P位点
	with open(plp, 'r') as f:
		for line in f:
			linenum += 1
			if linenum == 1:
				continue
			tmp = line.replace('\n', '').split('\t')
			site = tmp[0] + '_' + tmp[1] + '_' + tmp[2] + '_' + tmp[3]
			sites_p[site] = 'Pathogenic'

	## B/LB位点
	with open(blb, 'r') as f:
		for line in f:
			linenum += 1
			if linenum == 1:
				continue
			tmp = line.replace('\n', '').split('\t')
			site = tmp[1] + '_' + tmp[2] + '_' + tmp[3] + '_' + tmp[4]
			sites_b[site] = 'Benign'

	# print('sites length====', len(sites))
	sampleIDs = []
	nodeList = []
	totalList = []  # 以样本为单位，记录每个样本中所有变异的基因型
	varGTList = {}  # 以位点为单位，记录每个位点所有样本的基因型
	annoList = {} #记录每个位点的VEP注释信息
	with open(input_tsv, 'r') as f:  # open sourcefile
		lineNum = 0
		for line in f:
			lineNum += 1
			tmp = line.replace('\n', '').split('\t')
			# print(tmp)
			if (lineNum == 1):
				sampleIDs = list(tmp[12:])  # sampleIDs of GT
				# print(sampleIDs)
				for i in range(len(sampleIDs)):
					totalList.append([])  # define list to save GT
			if tmp[5] != feature:
				continue

			snp = tmp[0] + '_' + tmp[1] + '_' + tmp[2] + '_' + tmp[3]
			#记录位点的VEP信息，后续添加到网络中，作为节点的属性
			Symbol = tmp[4]
			Transcript = tmp[5]
			VEP_IMPACT = tmp[6]
			VEP_Consequence = tmp[8]
			AF_Max = tmp[9]
			VEP_HGVSc = tmp[10]
			VEP_HGVSp = tmp[11]


			# print(snp)
			label0 = 'V_'
			# for i in sites:
			#     if (i == snp):
			#         label0 = 'P_'
			for site, label in sites_p.items():
				if site == snp:
					if 'P' in label:  # 大写的P，以Pathogenic开头：Pathogenic, Pathogenic/Likely_pathogenic
						label0 = 'P_'
					else:  # Likely_pathogenic
						label0 = 'LP_'
				else:
					continue
			for site, label in sites_b.items():
				if site == snp:
					if 'Benign' in label:
						label0 = 'B_'
				else:
					continue
			vus1 = label0 + snp
			# print(tmp[0])
			sampleTmp = list(tmp[12:])  # get GT per line
			mod = 0
			for i in sampleTmp:
				if i == '1/1':
					mod += 1
			if(mod==1):
				label1 = '_0.5'
			elif(mod>1):
				label1 ='_1'
			else:
				label1 ='_0'

			nodeList.append(vus1+label1)

			annoList[vus1+label1] = {
				'Symbol': Symbol,
				'Transcript': Transcript,
				'VEP_IMPACT': VEP_IMPACT,
				'VEP_Consequence': VEP_Consequence,
				'AF_Max': AF_Max,
				'VEP_HGVSc': VEP_HGVSc,
				'VEP_HGVSp': VEP_HGVSp
			}
			# print(sampleTmp)
			for i in range(len(sampleTmp)):  # save GT to tatalList
				totalList[i].append(sampleTmp[i])
			varGTList[vus1+label1] = sampleTmp  # 位点在所有样本中的基因型
	# transList中的位点打上P/LP/V标签
	for i in range(len(transList)):
		var1 = list(transList[i])[0]
		var2 = list(transList[i])[1]
		if var1 in sites_p.keys():
			if 'P' in sites_p[var1]:
				var1 = 'P_' + var1
			else:
				var1 = 'LP_' + var1
		else:
			var1 = 'V_' + var1
		if var2 in sites_p.keys():
			if 'P' in sites_p[var1]:
				var2 = 'P_' + var2
			else:
				var2 = 'LP_' + var2
		else:
			var2 = 'V_' + var2
		transList[i] = frozenset([var1, var2])

	return nodeList, sampleIDs, totalList, varGTList, transList, annoList

def network(tsv, plp, blb,transcript, transList):
	# print('before',transList)
	nodeList, sampleIDs, totalList, varGTList, transList , annoList = labeling(tsv, plp, blb, transcript, transList)
	nodeNum = len(nodeList)
	print('nodeNum', nodeNum)
	# print('after',transList)
	filename = os.path.basename(tsv)
	smplst_x = filename.split('_')[-1].split('.')[0]

	# 对每个样本进行操作--以样本为单位，记录其发生<!杂合> 变异的位点
	sampleDict = {} #<!杂合>
	for i in range(len(totalList)):
		sampleID = sampleIDs[i]  # 当前样本的sampleID
		colList = totalList[i]  # 当前样本的所有位点基因型
		for j in range(nodeNum):  # 对每个位点进行操作
			nodeStr = nodeList[j]  # 当前位点的标记信息，如'V_13_20763104_T_C_0'
			if colList[j] == '0/1' or colList[
				j] == '1/0':  # 如果当前样本的当前位点为<!杂合>变异位点，则将该位点标记储存在nodeNumDict中：key为位点标记，value为1，形如{'P_13_20763553_CA_C_0': 1}
				# print(colList[j])
				if sampleID not in sampleDict.keys():  # 如果当前样本尚未记录在sampleDict字典中，则记录nodeNumDict到sampleDict字典中
					nodeNumDict = {}
					nodeNumDict[nodeStr] = 1
					sampleDict[sampleID] = nodeNumDict
				else:  # 如果当前样本已经在sampleDict字典中有记录，则新增该样本的位点标记
					nodeNumDict = sampleDict[sampleID]
					if nodeStr not in nodeNumDict.keys():  # 如果该位点尚未被记录，则添加当前位点
						nodeNumDict[nodeStr] = 1
					else:
						nodeNumDict[nodeStr] += 1  # ? 应该不会，一个样本的一个位点只会遍历一次。值不会是2？
					sampleDict[sampleID] = nodeNumDict
			else:
				continue
	# 对每个位点进行操作-以位点为单位，记录其发生纯合变异的样本; 对每个位点进行操作-以位点为单位，记录发生变异的位点

	varHomDict = {}  # 纯合
	varMutDict = [] #突变
	for i in range(len(varGTList)):
		varID = nodeList[i] #当前位点
		# print('varGTList[varID]===',varGTList[varID])
		gtList = varGTList[varID] #当前位点的所有样本的基因型
		for j in range(len(sampleIDs)): #对每个样本进行操作
			smpStr = sampleIDs[j] #当前样本ID
			if gtList[j] == '1/1':
				varMutDict.append(varID)
				if varID not in varHomDict.keys(): #当前位点尚未被记录，则添加记录
					smpDict = []
					smpDict.append(smpStr)
					varHomDict[varID] = smpDict
				else:  #当前位点已经在varHomDict字典中记录，则新增该位点的纯合突变样本ID
					if smpStr not in varHomDict[varID]: #当前样本尚未被记录，则添加该样本
						smpDict = varHomDict[varID]
						smpDict.append(smpStr)
						varHomDict[varID] = smpDict
			elif gtList[j] == '1/0' or gtList[j] == '0/1':
				varMutDict.append(varID)
			else:
				continue

	# 输出在单个样本中有同时发生<!杂合> 变异的位点
	# with open(f'{transcript}_multiVars.tsv', 'w') as f:
	# 	f.write('SampleID' + '\t' + 'CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'Transcript' + '\n')
	for sample, varsDict in sampleDict.items():
		vars = list(varsDict.keys())
		if len(vars) > 1:
			for i in range(len(vars)):
				var_tmp = sample + '\t' + '\t'.join(vars[i].split('_')[1:]) + '\t' + transcript
				# print('vars_tmp===',var_tmp)
				with open(f'{transcript}_multiVars.tsv', 'a') as f:
					f.write(var_tmp + '\n')

	# 输出与某位点同时发生突变的位点、样本，weight，phase,score

	varsFreqDict = {}  # 记录两个位点同时在一个样本中发生变异的频次信息，形如：{frozenset({'V_chr1_103005900_CCATCAT_CCAT', 'V_chr1_102914823_TAA_TAAA'}): 9
	for sample, varsDict in sampleDict.items():
		vars = list(varsDict.keys())
		print('vars==========', vars)
		smp = sample.split(':')[0].split(']')[1]  # 仅保留样本名，[26]AH-021:GT --> AH-021

		for i in range(len(vars)):
			node1 = vars[i]
			print('node1==========', node1)
			for j in range(i + 1, len(vars)):
				node2 = vars[j]
				print('node2==========', node2)
				if frozenset([node1, node2]) not in varsFreqDict.keys():
					freqSmpDict = []
					freqSmpDict.append(smp)
					varsFreqDict[frozenset([node1, node2])] = freqSmpDict
					print('if varsFreqDict[frozenset([node1,node2])]==========',
						  varsFreqDict[frozenset([node1, node2])])
				else:
					if smp not in varsFreqDict[frozenset([node1, node2])]:
						freqSmpDict = varsFreqDict[frozenset([node1, node2])]
						freqSmpDict.append(smp)
						varsFreqDict[frozenset([node1, node2])] = freqSmpDict
					print('else varsFreqDict[frozenset([node1,node2])]==========',
						  varsFreqDict[frozenset([node1, node2])])
				print('wai varsFreqDict[frozenset([node1,node2])]==========', varsFreqDict[frozenset([node1, node2])])

	print('varsFreqDict====', varsFreqDict)

	# 构建变异位点矩阵，遍历varHomDict和varsFreqDict，设矩阵值,设为频次值
	nodeNum = len(nodeList)
	# print('nodeNum',nodeNum)
	matrix = np.zeros((nodeNum, nodeNum), dtype=[('value', int), ('sample', object)])

	# 1)纯合变异位点——对角线
	for var, smps in varHomDict.items():
		index = nodeList.index(var)
		matrix[index, index] = (len(smps), smps)

	# 2)复合杂合变异位点——非对角线
	for varPair, smps in varsFreqDict.items():
		# print('varPair====',varPair)
		# print('freq====',freq)
		var1 = list(varPair)[0]
		var2 = list(varPair)[1]
		index1 = nodeList.index(var1)  # 查找node1的索引位置
		index2 = nodeList.index(var2)

		matrix[index1, index2] = (len(smps), smps)
		matrix[index2, index1] = (len(smps), smps)
		# print('matrix=====', matrix[index1,index2]['value'], matrix[index1,index2]['sample'])
		# print('type====',type( matrix[index1,index2]['sample']))

	# in trans文件中的位点是否记录在矩阵中，否，添加记录
	for i in range(len(transList)):
		var1 = list(transList[i])[0]
		var2 = list(transList[i])[1]
		if var1 in nodeList and var2 in nodeList:
			index1 = nodeList.index(var1)
			index2 = nodeList.index(var2)
			# print('before matrix[index1,index2]===',matrix[index1,index2])
			if matrix[index1, index2]['value'] == 0:
				matrix[index1, index2]['value'] = 1
				matrix[index2, index1]['value'] = 1
				# print('after matrix[index1,index2]===', matrix[index1, index2])

	G = nx.Graph()  # 构建网络
	#将位点*位点矩阵转换为图, 矩阵值不为0的点，均构建在图中。在同一个样本中有2个及以上<杂合>变异的位点,并用边连接，再加上权重; 纯合-自连接
	for i in range(nodeNum):
		# print('i=====', i)
		for j in range(nodeNum):
			# print('j======', j)
			# print('matrix[i][j] =======', matrix[i][j])
			if matrix[i][j]['value'] != 0:
				# print(i, j, matrix[i][j])
				# print('nodeList[i]======', nodeList[i])
				# print('nodeList[j]=========', nodeList[j])
				# 两个位点情况：P+V/P+P/V+V，根据P个数，设置边属性，后续画图时根据该属性值，决定是否画V+V
				label1_nodei = nodeList[i].split('_')[0]
				label1_nodej = nodeList[j].split('_')[0]
				pNode = 0
				if i != j : #复合杂合
					if label1_nodei == 'P':
						pNode += 1
					else:
						if label1_nodej == 'P':
							pNode += 1
				else: #纯合，自连接
					if label1_nodei == 'P':
						pNode += 1
				# G.add_edge(nodeList[i], nodeList[j], weight=matrix[i][j]['value'], sample=matrix[i][j]['sample'], pNode = pNode)
				#sample 列表转换为字符串
				sample_str = json.dumps(matrix[i][j]['sample']) if isinstance(matrix[i][j]['sample'], (list, dict)) else matrix[i][j]['sample']
				G.add_edge(nodeList[i], nodeList[j], weight=matrix[i][j]['value'], sample=sample_str,pNode=pNode)
				# 添加节点信息
				nx.set_node_attributes(G,{ nodeList[i] : annoList[nodeList[i]] })
				nx.set_node_attributes(G, { nodeList[j] : annoList[nodeList[j]] })
				# print('G=====', G)
	#未在图中的单个变异的位点加入到图中
	for i in varMutDict:
		if i not in G.nodes(): #未在图中
			G.add_node(i)
			nx.set_node_attributes(G,{ i : annoList[i] })

	wDegs = []
	Degs = []
	for node in G.nodes():
		wDegs.append(G.degree(node, weight='weight'))
		Degs.append(G.degree())
	# 设置图的布局
	# pos = nx.spring_layout(G)  # 使用 spring 布局
	#
	# print('plt=========')
	# # 获取边的权重
	# edge_labels = nx.get_edge_attributes(G, "weight")
	#
	# nx.draw(G, with_labels=True, node_color='lightblue', edge_color='gray', node_size=80, font_size=4)
	# nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=4, label_pos=1)
	# # 在图上标注权重
	# print('plot network=====')
	# plt.savefig('network_'+ transcript + '.pdf')
	#  获取当前transcript对应的genesymbol
	symbol = next((v['Symbol'] for v in annoList.values() if v['Transcript'] == transcript),None)

	print('save netdata=====')
	# nx.write_graphml(G,'network_' + transcript  + '_' + smplst_x  +  '.graphml')
	nx.write_graphml(G,'network_' + symbol + '_' + transcript  +  '.graphml')
	#保存前需将图中list,dict等复杂类型转换为字符串；读取时，需要还原q
	# for node, data in G.nodes(data=True):
	# 	for key, value in data.items():
	# 		if isinstance(value, (list,dict)):
	# 			data[key] = json.dumps(value)
	# nx.write_graphml(G, 'network_' + transcript + '.graphml')
	#读取时，需要还原
	# G = nx.read_graphml("graph.graphml")
	# for node, data in G.nodes(data=True):
	# 	for key, value in data.items():
	# 		try:
	# 			data[key] = json.loads(value)
	# 		except (json.JSONDecodeError, TypeError):
	# 			pass

	# import pickle
	# with open('netdata.pkl', 'wb') as f:
	# 	pickle.dump(G, f)



	# import igraph as ig
	#
	# G_ig = ig.Graph.from_networkx(G)
	# layout = G_ig.layout("fr")  # Fruchterman-Reingold布局
	# ig.plot(G_ig,
	#      layout=layout,
	#      vertex_color=G_ig.vs["color"],
	#      edge_width=G_ig.es["weight"],
	#      vertex_label=G_ig.vs["name"])
	# plt.show()

	return G, nodeList, varGTList, transList


def scoring(transcript, tsv, plp, blb, output, trans=None):
	print('rrrrrrr', transcript, tsv, plp, output, trans)
	print('trans in scoring =====', trans)
	transList = []
	if trans is not None:
		# 记录in trans 文件中的位点
		sitesList = {}
		with open(trans, 'r') as f:
			linenum = 0
			for line in f:
				linenum += 1
				if linenum == 1:
					continue
				tmp = line.replace('\n', '').split(
					'\t')
				sampleID = tmp[0]
				var_tmp = tmp[1] + '_' + tmp[2] + '_' + tmp[3] + '_' + tmp[4]
				if sampleID not in sitesList.keys():
					vars = []
					vars.append(var_tmp)
					sitesList[sampleID] = vars
				else:
					sitesList[sampleID].append(var_tmp)
			print('sitesList===', sitesList)

		for sample, vars in sitesList.items():
			transList.append(frozenset(vars))
		print('transList====', transList)
		# print('transsites===',transsites)
	print('transList===', transList)

	# with open(f'{transcript}_{output}', 'w') as f:
	# 	f.write('CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'Point_Hom' + '\t' +  'Symbol' + '\t' + 'Transcript' + '\t' + 'VEP_IMPACT' + '\t' + 'VEP_Consequence' + '\t' + \
	# 			'AF_Max' + '\t' + 'VEP_HGVSc' + '\t' + 'VEP_HGVSp' + '\t' + 'HomozygousNum' + '\t' + 'Score' + '\t' + 'PM3_Level' + '\n')
	# with open('freqSmps.txt', 'w') as f:
	#     f.write('Var1' + '\t' + 'Var2' + '\t' + 'Sample' + '\t' + 'Weight' + '\n')

	# with open(f'{transcript}_allFreqSmps.txt', 'w') as f:
	# 	f.write(
	# 		'PV' + '\t' + 'CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t'  + 'Point_Hom' + '\t' + 'Transcript' + '\t' + 'Var2' + '\t' + 'Sample' + '\t' + 'Weight' + '\n')
	G, nodeList, varGTList, transList = network(tsv, plp, blb, transcript, transList)
	transsites = [value for fs in transList for value in fs]  # in trans文件中的位点，打上P/LP/V标签，存为一个数组
	# 对在网络中的位点进行打分
	checkMatrix = []
	for node in G.nodes():
		# print('node===',node)
		nodeMatrix = []

		#节点属性
		Symbol = G.nodes[node]['Symbol']
		VEP_IMPACT = G.nodes[node]['VEP_IMPACT']
		VEP_Consequence = G.nodes[node]['VEP_Consequence']
		AF_Max = G.nodes[node]['AF_Max']
		VEP_HGVSc = G.nodes[node]['VEP_HGVSc']
		VEP_HGVSp = G.nodes[node]['VEP_HGVSp']

		nbs = G.adj[node]  # 访问node的邻居节点和边  #获得该变异位点的邻居节点及其权重

		# print('nbs====',nbs)
		var1 = node
		# 该位点纯合突变样本数
		homNum = varGTList[node].count('1/1')
		source = float(node.split('_')[-1]) #纯合突变分数
		if var1 in transsites:  # 当前位点在intrans文件
			# print('node1===',var1)
			score = 0
			pNum = 0
			lpNum = 0
			transPairNum_V = 0
			transPairNum_P = 0
			transPairNum_LP = 0
			wgt_PLP_trans = 0
			wgt_V_trans = 0
			wgt_P = 0
			wgt_LP = 0

			for nb, wgt in nbs.items():  # 遍历邻居节点，是否有邻居节点在intrans文件中
				# print('nb===',nb)
				print('nb==', nb)
				print('wgt==', wgt['weight'])
				var2 = nb
				if {var1,
					var2} in transList:  # 当前邻居节点与当前节点为in trans（不能单独判断是否在in trans文件,要与当前节点匹配判断） ------> confirmed in trans
					print('TTTTTT')
					nbLabel = nb.split('_')[0]
					if nbLabel == 'P':
						# pNum += 1
						transPairNum_P += 1
						wgt_PLP_trans += wgt['weight']
					elif nbLabel == 'LP':
						# lpNum += 1
						transPairNum_LP += 1
						wgt_PLP_trans += wgt['weight']
					else:
						transPairNum_V += 1
						wgt_V_trans += wgt['weight']
					print('wgt_PLP_trans==', wgt_PLP_trans)
					print('wgt_V_trans===', wgt_V_trans)
				else:  # 当前邻居节点与当前节点不是in trans，统计邻居节点P/LP个数。需根据邻居节点中P/LP位点的总个数来区分打分方式
					nbLabel = nb.split('_')[0]
					if nbLabel == 'P':
						pNum += 1
						wgt_P += wgt['weight']
					elif nbLabel == 'LP':
						lpNum += 1
						wgt_LP += wgt['weight']
				if (pNum + lpNum + transPairNum_P + transPairNum_LP) <= 1:
					if (
							transPairNum_P + transPairNum_LP + transPairNum_V) == 0:  # 邻居节点总PLP个数<=1，且无已知的in trans---> phase unknown
						score += 0.5 * wgt_P + 0.25 * wgt_LP
					else:  # 邻居节点总PLP个数<=1，但有已知的in trans:需加上该位点的纯合变异个数分值（0，0.5，1）；V-V的分值最多0.5
						score = source + 0.5 * wgt_P + 0.25 * wgt_LP + wgt_PLP_trans * 1 + (
							0.5 if wgt_V_trans * 0.25 >= 0.5 else wgt_V_trans * 0.25)
				else:  # 邻居节点PLP个数>1, ---->假定 in trans，须加上该位点纯合变异个数分值（0，0.5，1）
					score = source + (wgt_P + wgt_LP) / 2 + wgt_PLP_trans * 1 + (
						0.5 if wgt_V_trans * 0.25 >= 0.5 else wgt_V_trans * 0.25)
					print('score1=====', score)
				if score > 0 and score < 1:
					PM3_label = 'PM3_Supporting'
				elif score >= 1 and score < 2:
					PM3_label = 'PM3_Moderate'
				elif score >= 2 and score < 4:
					PM3_label = 'PM3_Strong'
				elif score >= 4:
					PM3_label = 'PM3_VeryStrong'
				else:
					PM3_label = 'NA'
				with open(f'{transcript}_{output}', 'a') as f:
					f.write('\t'.join(node.split('_')[1:]) + '\t' + Symbol + '\t' + transcript + '\t' + VEP_IMPACT + '\t' + VEP_Consequence + '\t' + \
							AF_Max + '\t' + VEP_HGVSc + '\t' + VEP_HGVSp + '\t' + str(homNum) + '\t' + str(score) + '\t' + str(PM3_label) + '\n')
		else:  # 当前位点不在in trans文件,判断邻居节点P个数
			pNum = 0
			lpNum = 0
			score = 0
			wgt_P = 0
			wgt_LP = 0
			print('node==', node)
			# alt = node.split('_')[-1]  # 一个变异单独输出一个文件，* 不能用于文件命名，需要特殊处理
			# if alt == '*':
				# cfn = node.split('_*')[0] + '_star'
			# else:
				# cfn = node
			# with open(cfn + '.txt', 'w') as f:
				# f.write('Var1' + '\t' + 'Var2' + '\t' + 'Sample' + '\t' + 'Weight' + '\n')
			for nb, wgt in nbs.items():  # 统计node的邻居节点nb和边wgt #初始化pNum=0，如果邻居节点为P标记的变异位点，则pNum+1
				print('nbs.items:', nb, wgt)

				# 输出与该位点同时发生变异的位点、样本信息
				# nodeMatrix.append(node + '\t' + nb + '\t' + ','.join(wgt['sample']) + '\t' + str(wgt['weight']) + '\t' + None + '\t' + 'None' +'\n')
				# with open(cfn + '.txt', 'a') as f:
					# f.write(node + '\t' + nb + '\t' + ','.join(wgt['sample']) + '\t' + str(wgt['weight']) + '\n')
					# f.write(node + '\t' + nb + '\t' + wgt['sample'] + '\t' + str(wgt['weight']) + '\n')
				with open(f'{transcript}_allFreqSmps.txt', 'a') as f:
					# f.write('\t'.join(node.split('_')) + '\t' + transcript + '\t' + nb + '\t' + ','.join(wgt['sample']) + '\t' + str(wgt['weight']) + '\n')
					f.write('\t'.join(node.split('_')) + '\t' + transcript + '\t' + nb + '\t' + wgt['sample'] + '\t' + str(wgt['weight']) + '\n')

				if node != nb:  # 非自连接，统计邻居节点P/LP个数
					nbLabel = nb.split('_')[0]  # 邻居节点的P/LP/V标记
					print('nbLabel', nbLabel)
					if nbLabel == 'P':
						pNum += 1
						wgt_P += wgt['weight']
					elif nbLabel == 'LP':
						lpNum += 1
						wgt_LP += wgt['weight']

				else: #自连接，跳过，避免将自身当成邻居节点，而错误计算，特别是单个自连接P。
					continue

				# print('type2222', type(','.join(wgt['sample'])), ','.join(wgt['sample']))
				# print('type333',wgt['weight'], type(wgt['weight']))
				# line = nb + '\t' + ','.join(wgt['sample']) + '\t' + str(wgt['weight'])+'\n'
				# print('line===',line,type(line))


				print('in for nbs.items:', pNum, lpNum, wgt_P, wgt_LP)
			print('over for nbs.items:', pNum, lpNum, wgt_P, wgt_LP)
			# 根据邻居节点P个数: 小于等于1时，Phase Unknown,根据邻居节点P/LP属性*权重统计得分+该位点纯合变异分数; 大于1时，假定in trans，1分*权重，再加上该位点的纯合变异个数得分值(0，0.5，1)
			if (pNum + lpNum) <= 1:
				score += 0.5 * wgt_P + 0.25 * wgt_LP + source
			else:  #权重不除以2
				score = (wgt_LP + wgt_P) + source

			# print('score2=====',score)
			print('score==', score)
			if score > 0 and score < 1:
				PM3_label = 'PM3_Supporting'
			elif score >= 1 and score < 2:
				PM3_label = 'PM3_Moderate'
			elif score >= 2 and score < 4:
				PM3_label = 'PM3_Strong'
			elif score >= 4:
				PM3_label = 'PM3_VeryStrong'
			else:
				PM3_label = 'NA'
			with open(f'{transcript}_{output}', 'a') as f:
				f.write('\t'.join(node.split('_')[1:]) + '\t' + Symbol + '\t' +  transcript + '\t' + VEP_IMPACT + '\t' + VEP_Consequence + '\t' + \
						AF_Max + '\t' + VEP_HGVSc + '\t' + VEP_HGVSp + '\t' + str(homNum)  + '\t' + str(score)+ '\t' + str(PM3_label) + '\n')

if __name__ == '__main__':

	print('sys.argv', len(sys.argv), sys.argv)
	print('sys.argv[1]', sys.argv[1])
	transcript = sys.argv[1]
	tsv = sys.argv[2]
	plp = sys.argv[3]
	blb = sys.argv[4]
	output = sys.argv[5]

	if len(sys.argv) == 6:
		intrans = None
	else:
		intrans = sys.argv[6]
	scoring(transcript, tsv, plp, blb, output, intrans)

	print('intrans', intrans)

