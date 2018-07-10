from rest_framework import permissions, generics
from rest_framework.decorators import api_view, permission_classes, detail_route, list_route
from rest_framework import mixins
from compounds import serializers
from compounds.models import *
from dynamic_rest import viewsets
from rest_framework.response import Response
from django_rdkit.models import *
from rest_framework.permissions import AllowAny
import sea
import pandas as pd
import smtplib
from email.mime.text import MIMEText
from email.header import Header
from smtplib import SMTP_SSL
from email.mime.multipart import MIMEMultipart
import os
# import pybel
from rdkit import Chem
import threading


TARGET_LIST = ['CHEMBL2034',
               'CHEMBL2035',
               'CHEMBL211',
               'CHEMBL4296',
               'CHEMBL275',
               'CHEMBL1833',
               'CHEMBL1836',
               'CHEMBL3729',
               'CHEMBL3535',
               'CHEMBL3594',
               'CHEMBL3746',
               'CHEMBL3119',
               'CHEMBL261',
               'CHEMBL3912',
               'CHEMBL4619',
               'CHEMBL1881',
               'CHEMBL4884',
               'CHEMBL252',
               'CHEMBL251',
               'CHEMBL256',
               'CHEMBL254',
               'CHEMBL1951',
               'CHEMBL4235',
               'CHEMBL3510',
               'CHEMBL4425',
               'CHEMBL245',
               'CHEMBL246',
               'CHEMBL241',
               'CHEMBL1940',
               'CHEMBL1942',
               'CHEMBL3977',
               'CHEMBL230',
               'CHEMBL2056',
               'CHEMBL3710',
               'CHEMBL1867',
               'CHEMBL3242',
               'CHEMBL3969',
               'CHEMBL4789',
               'CHEMBL232',
               'CHEMBL2973',
               'CHEMBL3805',
               'CHEMBL2609',
               'CHEMBL4409',
               'CHEMBL4408',
               'CHEMBL3012',
               'CHEMBL205',
               'CHEMBL229',
               'CHEMBL222',
               'CHEMBL223',
               'CHEMBL220',
               'CHEMBL221',
               'CHEMBL226',
               'CHEMBL224',
               'CHEMBL225',
               'CHEMBL2072',
               'CHEMBL1783',
               'CHEMBL4716',
               'CHEMBL1785',
               'CHEMBL2652',
               'CHEMBL291',
               'CHEMBL290',
               'CHEMBL1916',
               'CHEMBL3025',
               'CHEMBL1914',
               'CHEMBL3878',
               'CHEMBL2885',
               'CHEMBL217',
               'CHEMBL216',
               'CHEMBL214',
               'CHEMBL213',
               'CHEMBL4640',
               'CHEMBL5932',
               'CHEMBL210',
               'CHEMBL288',
               'CHEMBL4187',
               'CHEMBL286',
               'CHEMBL3231',
               'CHEMBL2326',
               'CHEMBL1900',
               'CHEMBL1821',
               'CHEMBL5163',
               'CHEMBL1827',
               'CHEMBL5451',
               'CHEMBL3421',
               'CHEMBL1987',
               'CHEMBL4267',
               'CHEMBL1980',
               'CHEMBL2717']

class DrugViewSet(viewsets.DynamicModelViewSet):
    queryset = Drug.objects.all()
    serializer_class = serializers.DrugSerializer

    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.data['smiles'])
        similarity = float(request.data['similarity'])
        substructure_search = int(request.data['substructure_search'])
        #perform substructure
        print(smiles, similarity, substructure_search)
        result = {}
        if substructure_search == 1:
            result = Drug.objects.filter(mol__hassubstruct=QMOL(Value(smiles))).all()
        # structure search
        else:
            try:
                result = Drug.objects.structure_search(smiles, similarity)
            except:
                print 'structure search error'

        if result:
            page = self.paginate_queryset(result)
            if page is not None:
                serializers = self.get_serializer(page, many=True)
                return self.get_paginated_response(serializers.data)
            serializer = self.get_serializer(result, many=True)
            return Response(serializer.data)
        return Response(result)

class DrugBankIDViewSet(viewsets.DynamicModelViewSet):
    queryset = DrugBankID.objects.all()
    serializer_class = serializers.DrugBankIDSerializer

    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.data['smiles'])
        similarity = float(request.data['similarity'])
        substructure_search = int(request.data['substructure_search'])
        # perform substructure
        print(smiles, similarity, substructure_search)
        result = {}
        if substructure_search == 1:
            result = DrugBankID.objects.filter(mol__hassubstruct=QMOL(Value(smiles))).all()
        # structure search
        else:
            try:
                result = DrugBankID.objects.structure_search(smiles, similarity)
            except:
                print 'structure search error'

        if result:
            page = self.paginate_queryset(result)
            if page is not None:
                serializers = self.get_serializer(page, many=True)
                return self.get_paginated_response(serializers.data)
            serializer = self.get_serializer(result, many=True)
            return Response(serializer.data)
        return Response(result)


class PathwayViewSet(viewsets.DynamicModelViewSet):
    queryset = Pathway.objects.all()
    serializer_class = serializers.PathwaySerializer

class TargetViewSet(viewsets.DynamicModelViewSet):
    queryset = Target.objects.all()
    serializer_class = serializers.TargetSerializer

class ChEMBL_small_molecule_all_infoViewSet(viewsets.DynamicModelViewSet):
    queryset = ChEMBL_small_molecule_all_info.objects.all()
    serializer_class = serializers.ChEMBL_small_molecule_all_info_Serializer

    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.data['smiles'])
        similarity = float(request.data['similarity'])
        substructure_search = int(request.data['substructure_search'])
        # perform substructure
        print smiles, similarity, substructure_search
        result = {}
        if substructure_search == 1:
            result = ChEMBL_small_molecule_all_info.objects.filter(mol__hassubstruct=QMOL(Value(smiles))).all()

        # structure search
        else:
            try:
                result = ChEMBL_small_molecule_all_info.objects.structure_search(smiles, similarity)
            except:
                print 'structure search error'

        if result:
            page = self.paginate_queryset(result)
            if page is not None:
                serializer = self.get_serializer(page, many=True)
                return self.get_paginated_response(serializer.data)
            serializer = self.get_serializer(result, many=True)
            return Response(serializer.data)

        return Response(result)


class ChEMBL_small_moleculeViewSet(viewsets.DynamicModelViewSet):
    queryset = ChEMBL_small_molecule.objects.all()
    serializer_class = serializers.ChEMBL_small_molecule_Serializer

    @list_route(methods=['POST', 'GET'], permission_classes=[permissions.AllowAny])
    def search(self, request):
        smiles = str(request.data['smiles'])
        similarity = float(request.data['similarity'])
        substructure_search = int(request.data['substructure_search'])
        # perform substructure
        print smiles, similarity, substructure_search
        result = {}
        if substructure_search == 1:
            result = ChEMBL_small_molecule.objects.filter(mol__hassubstruct=QMOL(Value(smiles))).all()
        # structure search
        else:
            try:
                result = ChEMBL_small_molecule.objects.structure_search(smiles, similarity)
            except:
                print 'structure search error'
        if result:
            page = self.paginate_queryset(result)
            if page is not None:
                serializer = self.get_serializer(page, many=True)
                return self.get_paginated_response(serializer.data)
            serializer = self.get_serializer(result, many=True)
            return Response(serializer.data)
        return Response(result)

@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def target_pred(request):
    smiles = str(request.data['smiles'])
    pred_data = sea.pred2_build_in_function(smiles)
    # pred_data = sea.pred2(smiles)
    print pred_data
    return Response(pred_data)

# @api_view(['POST'])
# @permission_classes([permissions.AllowAny])
# def bulk_target_pred(request):
#     structure_file = request.data['structure_file']
#     email_addr = str(request.data['email_addr'])
#     # f = open(structure_file)
#     lines = structure_file.readlines()
#     # f.close()
#     # print os.path.abspath(__file__)
#     # print os.path.dirname(os.path.abspath(__file__))
#     prediction_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pred_result', email_addr+'.csv')
#     df = pd.DataFrame()
#     # df.to_csv(prediction_result)
#     smile_list = []
#     pred_value_list = []
#     for idx, line in enumerate(lines):
#         line = line.strip()
#         pred_data = sea.pred2(line)
#         pred_data = [str(el) for el in pred_data]
#         smile_list.append(line)
#         pred_value_list.append(', '.join(pred_data))
#     df['smiles'] = smile_list
#     df['prediction'] = pred_value_list
#     df.to_csv(prediction_result)
#
#     host_server = 'smtp.qq.com'
#     sender_mail_addr = '1032847174@qq.com'
#     pwd = 'kyvixkkxcfadbccd'
#     receiver_mail_addr = email_addr
#     mail_content = 'Hello, this is the prediction result of your subscribed molecule'
#     mail_title = "JianpingLin's email"
#     msg = MIMEMultipart()
#     msg['Subject'] = Header(mail_title, 'utf-8')
#     msg['From'] = sender_mail_addr
#     msg['To'] = Header('Receiver', 'utf-8')
#
#     msg.attach(MIMEText(mail_content, 'html', 'utf-8'))
#     att1 = MIMEText(open(prediction_result).read(), 'base64', 'utf-8')
#     att1['Content-Type'] = 'application/octet-stream'
#     att1['Content-Disposition'] = 'attachment; filename="prediction_result"'
#     msg.attach(att1)
#
#     smtp = SMTP_SSL(host_server)
#     smtp.set_debuglevel(1)
#     smtp.ehlo(host_server)
#     smtp.login(sender_mail_addr, pwd)
#     smtp.sendmail(sender_mail_addr, receiver_mail_addr, msg.as_string())
#     smtp.quit()
#     return Response({})

def send_prediction_result(structure_file, email_addr, prediction_result, prediction_type):
    """
    :param prediction_type: users want to predict all target protein or one single.
    :return:
    """
    df = pd.DataFrame()
    if structure_file.name.endswith('sdf'):
        mols = Chem.ForwardSDMolSupplier(structure_file)
        for idx, mol in enumerate(mols):
            # print(mol.GetNumAtoms())
            if idx < 100:
                if mol is None:
                    continue
                smile = Chem.MolToSmiles(mol)
                if prediction_type == 'All':
                    pred_data = sea.pred2_all_target_list(smile)
                else:
                    pred_data = sea.pred2_single_target(smile, target=prediction_type)
                temp_df = pd.DataFrame(pred_data)
                temp_df['index'] = [idx] * len(pred_data)
                df = pd.concat([df, temp_df])
            else:
                break
    if structure_file.name.endswith('smi'):
        lines = structure_file.readlines()
        for idx, line in enumerate(lines):
            if idx < 100:
                line = line.strip()
                if prediction_type == 'All':
                    pred_data = sea.pred2_all_target_list(line)
                else:
                    pred_data = sea.pred2_single_target(line, target=prediction_type)
                temp_df = pd.DataFrame(pred_data)
                temp_df['index'] = [idx] * len(pred_data)
                df = pd.concat([df, temp_df])
            else:
                break

    # if structure_file.name.endswith('mol2'):
    #     # mols = Chem.MolFromMol2File(structure_file)
    #     str = structure_file.read()
    #     s = Chem.MolFromMol2Block()
    #     for mol in mols:
    #         print(mol.GetNumAtoms())
    df.fillna('None', inplace=True)
    print(df.columns)
    # df = df.reindex_axis(['index', 'smile'] + list(df.columns[:])-['index', 'smile'], axis=1)
    df = df[['index', 'smiles', 'atompair_hashed', 'maccs', 'morgan_hashed', 'topological_hashed', 'chembl_id']]
    df.to_csv(prediction_result)

    host_server = 'smtp.qq.com'
    sender_mail_addr = '1032847174@qq.com'
    pwd = 'kyvixkkxcfadbccd'
    receiver_mail_addr = email_addr
    mail_content = 'Hello, this is the prediction result of your subscribed molecule'
    mail_title = "JianpingLin's email"
    msg = MIMEMultipart()
    msg['Subject'] = Header(mail_title, 'utf-8')
    msg['From'] = sender_mail_addr
    msg['To'] = Header('Receiver', 'utf-8')

    msg.attach(MIMEText(mail_content, 'html', 'utf-8'))
    att1 = MIMEText(open(prediction_result).read(), 'base64', 'utf-8')
    att1['Content-Type'] = 'application/octet-stream'
    att1['Content-Disposition'] = 'attachment; filename="prediction_result.csv"'
    msg.attach(att1)

    smtp = SMTP_SSL(host_server)
    smtp.set_debuglevel(1)
    smtp.ehlo(host_server)
    smtp.login(sender_mail_addr, pwd)
    smtp.sendmail(sender_mail_addr, receiver_mail_addr, msg.as_string())
    smtp.quit()


@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def bulk_target_pred(request):
    structure_file = request.data['structure_file']
    email_addr = str(request.data['email_addr'])
    prediction_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pred_result', email_addr + '.csv')
    t = threading.Thread(target=send_prediction_result, args=(structure_file, email_addr, prediction_result))
    t.setDaemon(False)
    t.start()
    return Response('Congratulations, you have submitted the file successfully')


@api_view(['POST'])
@permission_classes([permissions.AllowAny])
def bulk_target_pred_2(request):

    """
    prediction_type: user can choonse all target protein or one single target protein (ChEMBL ID)
    """
    structure_file = request.data['structure_file']
    email_addr = str(request.data['email_addr'])
    prediction_type = str(request.data['prediction_type'])
    prediction_result = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pred_result', email_addr + '.csv')
    # send_prediction_result(structure_file, email_addr, prediction_result, prediction_type)
    t = threading.Thread(target=send_prediction_result, args=(structure_file, email_addr, prediction_result, prediction_type))
    t.setDaemon(False)
    t.start()
    return Response('Congratulations, you have submitted the file successfully. The prediction result will be sent to your email.')

class FeedbackCreateView(generics.CreateAPIView):
    queryset = Feedback.objects.all()
    serializer_class = serializers.FeedbackSerializer
    permission_classes = (AllowAny, )