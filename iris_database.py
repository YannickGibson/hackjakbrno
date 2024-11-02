import iris
import time
import os
import pandas as pd

class IrisDatabaseHandler:
    schema_name = 'DNASchema'
    vector_length = 768
    connection = None
    cursor = None

    def __init__(self):
        username = 'demo'
        password = 'demo'
        hostname = os.getenv('IRIS_HOSTNAME', 'localhost')
        port = '1972'
        namespace = 'USER'
        CONNECTION_STRING = f"{hostname}:{port}/{namespace}"

        #try except or preferably context manager
        self.connection = iris.connect(CONNECTION_STRING, username, password)
        self.cursor = self.connection.cursor()

    def check_table_exists(self, table_name: str) -> bool:
        self.cursor.execute(
            f"SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = '{self.schema_name}' AND table_name = '{table_name}'")

        return self.cursor.fetchone()[0]

    def insert_embedded_sequence(self, df: pd.DataFrame) -> None:
        table_name = f"Sequences"
        schema_table_name = f"{self.schema_name}.{table_name}"

        # print(str(df['sequence_vector'][0]))

        sql = f"""
            INSERT INTO {schema_table_name}
            (file_path, index, sequence_string, sequence_vector) 
            VALUES (?, ?, ?, TO_VECTOR(?, DOUBLE, {self.vector_length}))
            """

        start_time = time.time()
        # Prepare the list of tuples (parameters for each row)
        data = [
            (
                row['file_path'],
                row['index'],
                row['sequence_string'],
                str(row['sequence_vector'])
            )
            for index, row in df.iterrows()
        ]

        self.cursor.executemany(sql, data)
        end_time = time.time()
        print(f"time taken to add {len(df)} entries: {end_time - start_time} seconds")

    def search(self):
        table_name = f"Sequences"
        schema_table_name = f"{self.schema_name}.{table_name}"
        searchVector = "[0.005567946005612612, 0.05584058165550232, 0.05253725126385689, -0.1071929857134819, -0.09200528264045715, -0.08021321892738342, 0.018258729949593544, -0.053453296422958374, -0.06098167225718498, -0.0131150484085083, -0.10609208047389984, 0.10756513476371765, -0.04655534401535988, -0.1245153546333313, -0.021320588886737823, -0.024093730375170708, -0.011836654506623745, -0.0012545118806883693, -0.06228006258606911, -0.08083449304103851, 0.0838717371225357, 0.10319287329912186, -0.047839246690273285, 0.020847899839282036, 0.019623126834630966, 0.09284456074237823, 0.14692053198814392, 0.07748342305421829, -0.11310255527496338, -0.05824828892946243, -0.04241526871919632, 0.1120801642537117, -0.08281367272138596, -0.06172984838485718, -0.06301244348287582, 0.00458114966750145, -0.09086631238460541, -0.11640608310699463, 0.044896673411130905, -0.0022167591378092766, -0.07613535970449448, 0.059347912669181824, 0.06506941467523575, -0.02570527233183384, -0.01449514552950859, 0.18216344714164734, -0.19121041893959045, -0.061893630772829056, 0.00036359476507641375, 0.16914021968841553, 0.02291102707386017, -0.11774685233831406, -0.061718568205833435, 0.0711820125579834, -0.07138983905315399, 0.14886851608753204, 0.05848235636949539, -0.0831935927271843, -0.10326602309942245, 0.06968813389539719, -0.026936182752251625, 0.03825229033827782, -0.01862645521759987, -0.11848265677690506, 0.154217928647995, 0.06256397813558578, -0.07681677490472794, -0.05906852334737778, -0.01472118217498064, 0.12516556680202484, 0.10214132070541382, 0.08787070959806442, -0.06934914737939835, -0.10336168110370636, 0.2124495953321457, -0.07920413464307785, 0.023548630997538567, -0.045631930232048035, 0.06890986859798431, 0.0008616334525868297, -0.14392895996570587, 0.1281549334526062, 0.10233938694000244, 0.013151460327208042, -0.0015386701561510563, 0.08760885894298553, -0.025436945259571075, -0.12774312496185303, 0.1170009896159172, 0.051897577941417694, -0.026416223496198654, 0.027427537366747856, -0.1555347889661789, -0.08685432374477386, -0.06393925845623016, 0.11826861649751663, -0.037755727767944336, -0.0014531121123582125, 0.08900272846221924, -0.08142562955617905, -0.016928838565945625, 0.08317048102617264, -0.14606204628944397, -0.052445005625486374, 0.07228881865739822, -0.0034573369193822145, -0.15554703772068024, 0.060307830572128296, 0.08283442258834839, 0.03186560794711113, 0.10909057408571243, -0.05830824002623558, -0.09943757206201553, 0.028304964303970337, -0.07052242755889893, -0.03362653777003288, 0.020128363743424416, 0.05721117928624153, -0.01876099966466427, 0.09216048568487167, -0.1127840131521225, 0.08234231173992157, 0.002093872521072626, -0.07415980100631714, 0.03354586288332939, -0.2293737530708313, 0.16525429487228394, 0.023525122553110123, 0.11316780000925064, 0.004120740108191967, 0.07572415471076965, 0.07780192047357559, -0.01968574710190296, 0.04744898900389671, -0.013829611241817474, -0.07678254693746567, -0.016072377562522888, -0.03684350848197937, -0.17853327095508575, 0.06358842551708221, -0.08059680461883545, -0.07251589745283127, -0.027748797088861465, 0.08629864454269409, 0.045292120426893234, -0.00820470042526722, -0.12889353930950165, 0.0036877994425594807, -0.06192272529006004, 0.11671976000070572, -0.09121920168399811, -0.03675450012087822, 0.05176088958978653, 0.008361644111573696, 0.0013955743052065372, 0.05643901228904724, 0.07652964442968369, 0.05296171456575394, -0.07503605633974075, 0.11066874861717224, 0.08196873217821121, 0.02629515342414379, 0.036117713898420334, -0.035912033170461655, 0.06957980245351791, -0.09213823080062866, 0.11518070846796036, 0.09613464772701263, -0.0674707442522049, -0.05382147058844566, -0.13069720566272736, 0.05205922573804855, -0.060412563383579254, 0.1329909861087799, -0.11985227465629578, 0.09950721263885498, 0.0684075579047203, -0.1369023472070694, -0.10107140243053436, -0.024485060945153236, 0.15839943289756775, -0.06938152760267258, 0.13705024123191833, -0.1787838488817215, 0.10876212269067764, 0.06222466379404068, 0.07887239754199982, -0.04152544215321541, -0.17529776692390442, -0.02673589065670967, 0.04733636602759361, 0.12834972143173218, 0.036938123404979706, 0.0956173986196518, -0.0860091894865036, -0.15586940944194794, 0.000443686090875417, 0.026570983231067657, 0.02523086592555046, -0.08833930641412735, 0.09651626646518707, 0.027832835912704468, 0.08439339697360992, 0.10263197124004364, 0.1634080410003662, -0.11971070617437363, 0.039109714329242706, 0.04844221472740173, -0.09199974685907364, -0.06563159823417664, 0.08719392865896225, -0.016355769708752632, -0.2531794011592865, 0.016486050561070442, 0.03487856686115265, 0.12485796958208084, 0.023848390206694603, -0.1078108474612236, 0.006283764727413654, -0.039964210242033005, 0.10403931885957718, 0.08623172342777252, -0.018556004390120506, -0.018429715186357498, -0.03621677681803703, -0.06735444068908691, -0.0713944211602211, 0.20688709616661072, -0.03641575574874878, 0.0034590887371450663, 0.1697010099887848, -0.06836260110139847, -0.06474994868040085, 0.2394859790802002, 0.00387177593074739, -0.015456345863640308, -0.07000859826803207, 0.11943824589252472, -0.030189743265509605, -0.08995173871517181, 0.1168728917837143, -0.04408863186836243, 0.019475286826491356, 0.09197891503572464, 0.023390885442495346, 0.03764539584517479, -0.03412734344601631, -0.023574380204081535, 0.015638619661331177, -0.031023042276501656, -0.1295434832572937, 0.013699392788112164, 0.047505129128694534, 0.04509749263525009, 0.12153427302837372, -0.021549934521317482, 0.029106277972459793, -0.13642427325248718, 0.0971067026257515, 0.011164594441652298, 0.07037096470594406, -0.07670522481203079, 0.005124831572175026, 0.17342284321784973, 0.0409855954349041, 0.006655171047896147, -0.08789469301700592, 0.00011499186075525358, 0.0048828888684511185, 0.1039244532585144, -0.0958014652132988, 0.0637771487236023, 0.07804162055253983, 0.012114152312278748, 0.017724448814988136, 0.1319924145936966, -0.04436670243740082, 0.04501376301050186, 0.06208150461316109, 0.021960662677884102, 0.08457939326763153, -0.08037029951810837, -0.11056569963693619, 0.16797000169754028, 0.060173943638801575, -0.15908093750476837, -0.010669907554984093, 0.19360271096229553, 0.08529982715845108, 0.04144861176609993, -0.01702849194407463, 0.12988050282001495, -0.043098460882902145, 0.026335882022976875, 0.10952868312597275, -0.03598373010754585, 0.033814191818237305, 0.006414839066565037, 0.07184112071990967, 0.07479589432477951, 0.058214593678712845, 0.022850625216960907, 0.08926180005073547, 0.02479800209403038, 0.08465404808521271, -0.10134384036064148, -0.01732776314020157, 0.016589634120464325, -0.009917241521179676, -0.018154576420783997, 0.013674499467015266, -0.009872830472886562, -0.03190907463431358, 0.0010286930482834578, -0.005017325282096863, 0.07084416598081589, 0.1470610499382019, 0.09800207614898682, 0.1329558789730072, -0.022930631414055824, 0.06093785539269447, 0.11566345393657684, -0.03325872868299484, -0.039801593869924545, -0.07711540907621384, 0.022326502948999405, 0.16072003543376923, 0.06153549998998642, 0.043550241738557816, -0.13552631437778473, 0.08121009916067123, 0.06468807905912399, -0.03103460557758808, -0.17129109799861908, -0.015488930977880955, 0.09889931231737137, 0.17789460718631744, -0.06983828544616699, -0.09689456224441528, -0.013047678396105766, -0.003974795341491699, -0.22226214408874512, 0.11173264682292938, 0.10448063910007477, 0.1090790405869484, 0.0763501450419426, -0.1366875320672989, -0.042554568499326706, -0.14332284033298492, -0.0888737365603447, 0.03733329474925995, -0.1397484391927719, -0.05925695225596428, -0.004829659126698971, 0.11102909594774246, 0.07300281524658203, 0.08138810098171234, -0.35261690616607666, -0.043211616575717926, 0.05141213908791542, 0.058312494307756424, -0.0315588116645813, -0.09786821156740189, -0.13346518576145172, 0.025366803631186485, 0.060363661497831345, -0.20385028421878815, -0.07631443440914154, -0.04615223407745361, 0.10724145174026489, -0.07829232513904572, 0.14194537699222565, -0.057909127324819565, -0.005725025665014982, -0.03476885333657265, 0.06615164130926132, -0.06258561462163925, 0.013864927925169468, -0.07631272077560425, -0.004434129688888788, 0.05333052575588226, 0.0628020167350769, 0.04807561635971069, -0.05616830661892891, 0.08465202897787094, 0.1456976979970932, -0.06108113378286362, 0.10006009042263031, -0.008152663707733154, 0.030165813863277435, 0.054266929626464844, 0.19770340621471405, -0.005060201045125723, 0.10113418102264404, -0.001262869336642325, -0.07456525415182114, -0.1772809475660324, 0.13545192778110504, 0.025097129866480827, 0.06746228039264679, -0.130597785115242, -0.045684002339839935, 0.026219310238957405, 0.018691498786211014, 0.12212701886892319, -0.03277352452278137, 0.084786057472229, 0.2365441471338272, 0.041564252227544785, 0.016925128176808357, 0.17611807584762573, 0.05103185772895813, 0.0074379718862473965, 0.022543281316757202, -0.06975095719099045, -0.11128799617290497, -0.08619343489408493, 0.09939354658126831, 0.16960100829601288, -0.09894761443138123, 0.12230093777179718, 0.08435428142547607, -0.11758396774530411, -0.09416625648736954, -0.006350617855787277, 0.06716731190681458, 0.035510897636413574, -0.061804160475730896, 0.04593326896429062, 0.0012577985180541873, -0.04252268746495247, 0.05302158743143082, 0.008476796559989452, 0.13365671038627625, -0.03842144459486008, -0.09872905910015106, 0.07861033082008362, 0.014285103417932987, -0.04748281091451645, -0.699170708656311, -0.0015231125289574265, 0.12593646347522736, 0.04867700859904289, -0.08005878329277039, 0.0652301236987114, 0.009427541866898537, 0.04631786793470383, -0.03986922651529312, -0.1125975027680397, -0.025472048670053482, 0.07853347063064575, -0.036765944212675095, 0.1310224086046219, 0.015152249485254288, -0.09234835207462311, 0.07913850247859955, 0.14671993255615234, -0.037818554788827896, 0.015937253832817078, -0.012735537253320217, 0.013893388211727142, -0.011100499890744686, -0.14689543843269348, -0.01390242762863636, -0.056651849299669266, 0.11532227694988251, -0.15534332394599915, -0.052005834877491, 0.07510270178318024, 0.09272515773773193, -0.060546875, -0.11128555983304977, 0.015481244772672653, 0.10663043707609177, -0.010663741268217564, -0.08495376259088516, -0.10890527814626694, 0.039998654276132584, 0.09846397489309311, -0.004053547512739897, -0.08005595207214355, -0.1454671174287796, 0.008917181752622128, 0.06333544850349426, 0.02473784238100052, -0.08282723277807236, 0.04385726526379585, 0.036424390971660614, 0.17036214470863342, 0.037702277302742004, 0.22596698999404907, 0.01681077852845192, -0.008757512085139751, -0.014721989631652832, 0.11991573125123978, 0.024072417989373207, 0.023691583424806595, 0.0644305944442749, 0.03182271122932434, 0.04265553876757622, 0.19966372847557068, 0.02056116610765457, -0.07269701361656189, 0.02482989989221096, 0.0741850957274437, 0.08220775425434113, -0.05475340038537979, -0.10792157053947449, -0.08832444250583649, -0.012940528802573681, -5.692243576049805e-06, -0.08539236336946487, 0.0034131044521927834, -7.380203169304878e-05, 0.13069291412830353, -0.07501109689474106, -0.03616290166974068, 0.1144862174987793, 0.08882800489664078, -0.07113967835903168, 0.1178596094250679, 0.0599922239780426, -0.01898936554789543, 0.11543125659227371, -0.10705562680959702, -0.16161879897117615, 0.1148393452167511, 0.16864994168281555, 0.03788064047694206, 0.006618945859372616, -0.0978311076760292, 0.0791579857468605, 0.11583047360181808, 0.1082996055483818, -0.055451419204473495, -0.034564074128866196, -0.13431672751903534, -0.01981782168149948, -0.018322480842471123, 0.023129448294639587, -0.052597321569919586, -0.1014292761683464, -0.03806028515100479, -0.03278294950723648, 0.07066873461008072, -0.1322256475687027, -0.012916470877826214, 0.08478444069623947, -0.033756937831640244, -0.15508638322353363, -0.2935594916343689, 0.029988640919327736, 0.0387076660990715, 0.059736356139183044, -0.09624260663986206, 0.09071588516235352, 0.10449741780757904, 0.11523397266864777, 0.04039856418967247, 0.004885642323642969, -0.27447277307510376, 0.07372815161943436, 0.1162089928984642, 0.040372442454099655, 0.0110799390822649, 0.14609789848327637, -0.00874114315956831, 0.06324312835931778, -0.018592264503240585, 0.07578467577695847, 0.08644925057888031, -0.010339930653572083, -0.012324915267527103, -0.12391521781682968, -0.007716946303844452, -0.08021609485149384, -0.09614095091819763, 0.08081022650003433, 0.02809085138142109, -0.08864371478557587, -0.09327809512615204, -0.0249694362282753, 0.05840325355529785, 0.006842431146651506, 0.15638212859630585, 0.08296877890825272, -0.09303031116724014, -0.008194287307560444, 0.03563850000500679, -0.0727328211069107, 0.09536231309175491, -0.09367617964744568, 0.11526406556367874, 0.015146215446293354, 0.1221894845366478, 0.00376216066069901, 0.0449809655547142, 0.09813596308231354, 0.008344391360878944, 0.08394187688827515, -0.0656084194779396, 0.04916887730360031, 0.08341381698846817, -0.17892839014530182, -0.005237270146608353, -0.09894675016403198, 0.15810993313789368, -0.043968748301267624, 0.045117463916540146, -0.02105139195919037, 0.017994705587625504, -0.020199505612254143, -0.03385104238986969, 0.05343829467892647, 0.11353370547294617, -0.02525232546031475, -0.04555286094546318, -0.11228816211223602, 0.11476180702447891, 0.01842418685555458, -0.06678750365972519, 0.034736767411231995, -0.16359545290470123, 0.001152417971752584, 0.08800476044416428, 0.06312008202075958, -0.02025621570646763, 0.16132555902004242, -0.07756703346967697, 0.03127997741103172, 0.07318521291017532, -0.00999964214861393, -0.01672602817416191, 0.18560412526130676, 0.1093418300151825, 0.07500232756137848, -0.08282678574323654, 0.04810117185115814, -0.16061076521873474, -0.0758066475391388, -0.06679459661245346, -0.2835156321525574, -0.05939658358693123, -0.08143658936023712, -0.07754328101873398, 0.1092558279633522, 0.07895197719335556, -0.10230867564678192, -0.11468994617462158, -0.00633363937959075, -0.0715075209736824, -0.007727963849902153, -0.05273440480232239, -0.03265032544732094, -0.047487929463386536, -0.024892576038837433, -0.010919672437012196, 0.017844244837760925, -0.0922100767493248, -0.027980420738458633, -0.045575689524412155, -0.11165206134319305, -0.13106027245521545, 0.0637231096625328, 0.14181552827358246, 0.0054133799858391285, 0.02492358349263668, 0.09700897336006165, 0.045387350022792816, 0.1000889465212822, -0.06783052533864975, -0.05302988365292549, 0.07885540276765823, 0.08087228238582611, -0.04427128657698631, 0.03366127237677574, -0.03288160637021065, 0.006763693410903215, 0.033283691853284836, 0.07099156081676483, 0.11627204716205597, -0.1694198101758957, -0.0627366453409195, -0.172700434923172, 0.09012974798679352, -0.03269866108894348, 0.03003133274614811, -0.11981171369552612, -0.062138885259628296, -0.07539636641740799, 0.028097694739699364, -0.008289037272334099, 0.008683780208230019, 0.029693542048335075, -0.02305370382964611, -0.05663926526904106, -0.024823367595672607, 0.03923086076974869, -0.03885086625814438, 0.008908485993742943, -0.06741735339164734, -0.0786290392279625, -0.006367533467710018, 0.07953144609928131, 0.059036970138549805, -0.11863056570291519, 0.0394853875041008, 0.13067473471164703, -0.15410947799682617, -0.06403090804815292, -0.05203595757484436, 0.12690871953964233, 0.1190241202712059, 0.15231063961982727, -0.02086498960852623, -0.11984669417142868, 0.17170995473861694, -0.029949063435196877, -0.09570105373859406, 0.11324670910835266, 0.05100725218653679, 0.1095539778470993, -0.011797498911619186, -0.04561547562479973, 0.03441830724477768, 0.08028580993413925, -0.16597607731819153, 0.04089171811938286, -0.0831545814871788, -0.017229489982128143, 0.11600635200738907, -0.026016363874077797, -0.16478702425956726, 0.08125145733356476, 0.01585453189909458, -0.007084742654114962, 0.08825679868459702, -0.06859417259693146, -0.0806843638420105, -0.015784189105033875, 0.0902893990278244, 0.023817718029022217, -0.008258815854787827, -0.051332246512174606, -0.062020134180784225, 0.10455089062452316, 0.020050041377544403, 0.11334715783596039, -0.06829625368118286, 0.007472225930541754, 0.016630088910460472, -0.1491904854774475, 0.0706552043557167, 0.08364774286746979, 0.05164001137018204, 0.19535768032073975, -0.005448841955512762, 0.02117820829153061, 0.10171586275100708, -0.12072812020778656, -0.07307090610265732, 0.00031729767215438187, -0.1357852667570114, -0.04002280533313751, 0.009914297610521317, 0.051014676690101624, 0.10548384487628937, 0.04187123849987984]"

        sql = f"""
            SELECT TOP ? file_path, index, sequence_string, sequence_vector
            FROM {schema_table_name}
            ORDER BY VECTOR_DOT_PRODUCT(sequence_vector, TO_VECTOR(?)) DESC
        """

        numberOfResults = 3

        # Execute the query with the number of results and search vector as parameters
        self.cursor.execute(sql, [numberOfResults, str(searchVector)])

        # Fetch all results
        results = self.cursor.fetchall()
        for row in results:
            print(row)

    def create_tables(self):
        table_name = f"Sequences"
        schema_table_name = f"{self.schema_name}.{table_name}"
        # Possibly longtext for bigger strings sequences
        table_definition = f"(file_path VARCHAR(255), index INT, sequence_string TEXT, sequence_vector VECTOR(DOUBLE, {self.vector_length}))"


        if self.check_table_exists(table_name):
            print("Table already exists! Aborting...")
            return

        self.cursor.execute(f"CREATE TABLE {schema_table_name} {table_definition}")

    def drop_table(self):
        table_name = f"Sequences"
        schema_table_name = f"{self.schema_name}.{table_name}"

        if not self.check_table_exists(table_name):
            print("Table does not exist! Aborting...")
            return

        self.cursor.execute(f"DROP TABLE {schema_table_name}")

db = IrisDatabaseHandler()

# db.drop_table()
db.create_tables()
# db.search()

df = pd.read_csv('example_rows_for_db.csv', sep=';')
# print(df.head())
#
IrisDatabaseHandler().insert_embedded_sequence(df)
