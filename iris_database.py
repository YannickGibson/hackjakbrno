import iris
import time
import os
import pandas as pd
from tqdm import tqdm

class IrisDatabaseHandler:
    schema_name = 'DNASchema'
    vector_length = 768
    connection = None
    cursor = None
    results_number = 1

    table_name = f"Sequences"
    schema_table_name = f"{schema_name}.{table_name}"

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
        sql = f"""
            INSERT INTO {self.schema_table_name}
            (barcode, file_path, index, sequence_string, sequence_vector)
            VALUES (?, ?, ?, ?, TO_VECTOR(?, DOUBLE, {self.vector_length}))
            """

        start_time = time.time()
        data = [
            (
                row['barcode'],
                row['file_path'],
                row['index'],
                row['sequence_string'],
                str(row['sequence_vector'])
            )
            for index, row in tqdm(df.iterrows(), total=4000*9)
        ]

        self.cursor.executemany(sql, data)
        end_time = time.time()
        print(f"time taken to add {len(df)} entries: {end_time - start_time} seconds")

    def search(self, barcode: str, search_vector: str) -> list[tuple]:
        sql = f"""
            SELECT TOP ? barcode, file_path, index, sequence_string, sequence_vector
            FROM {self.schema_table_name}
            WHERE barcode=?
            ORDER BY VECTOR_COSINE(sequence_vector, TO_VECTOR(?)) DESC
        """

        # Execute the query with the number of results and search vector as parameters
        self.cursor.execute(sql, [self.results_number, barcode, search_vector])

        # Fetch all results
        return self.cursor.fetchall()

    def create_table(self):
        # Possibly longtext for bigger strings sequences
        table_definition = f"(barcode VARCHAR(255), file_path VARCHAR(255), index INT, sequence_string TEXT, sequence_vector VECTOR(DOUBLE, {self.vector_length}))"


        if self.check_table_exists(self.table_name):
            print("Table already exists! Aborting...")
            return

        self.cursor.execute(f"CREATE TABLE {self.schema_table_name} {table_definition}")

    def drop_table(self):
        if not self.check_table_exists(self.table_name):
            print("Table does not exist! Aborting...")
            return

        self.cursor.execute(f"DROP TABLE {self.schema_table_name}")


if __name__ == "__main__":
    db = IrisDatabaseHandler()
    #
    db.drop_table()
    db.create_table()
    # res = db.search("barcode20", "-.050322711467742919921, .10289134830236434936, .083920463919639587402, -.12247113138437271118, -.10290934145450592041, -.049007054418325424194, -.069650150835514068603, -.17487174272537231446, .034322228282690048217, .045333921909332275391, -.11466582864522933959, .042310345917940139771, .031600851565599441528, -.079703688621520996093, .083531059324741363526, -.021811420097947120666, -.057629127055406570434, .10324039310216903686, -.0034290526527911424636, -.12344066053628921508, .045306257903575897216, .12903387844562530517, .036231078207492828369, .091477788984775543212, .0072188684716820716857, .026461940258741378784, .29882004857063293457, .037163719534873962402, -.099755316972732543946, -.044139072299003601074, -.059011001139879226684, .083978697657585144042, -.034565322101116180419, .031997758895158767701, -.073293648660182952881, -.015557189472019672393, -.079128697514533996582, .043456509709358215332, -.034648526459932327271, .034782122820615768432, -.035960726439952850341, .032327037304639816284, -.076798893511295318603, .00075466901762410998344, -.029795778915286064147, .037644252181053161621, -.094308562576770782471, -.024658143520355224609, -.058196015655994415283, .072090260684490203857, .040109198540449142456, -.024893401190638542176, -.059397257864475250244, .031919520348310470581, -.057720433920621871948, .022258903831243515014, .11900207400321960449, -.17419515550136566162, -.12056382745504379272, -.028444508090615272521, .0051667802035808563232, .077476315200328826904, -.036257214844226837158, -.010945518501102924346, .10691558569669723511, .14685280621051788331, -.031189424917101860046, -.043201219290494918823, -.078069634735584259033, .022231418639421463012, .063057556748390197753, .12307672202587127686, -.037790060043334960937, -.063286721706390380859, .20479558408260345458, -.061841845512390136718, .0039520636200904846191, -.026175521314144134521, .063074663281440734863, -.063928104937076568603, -.031487476080656051636, -.014378922991454601287, .095356501638889312744, .011652008630335330963, .038628049194812774658, .049804642796516418457, .0044902404770255088806, -.12368037551641464233, .11772273480892181396, .037090521305799484252, -.17080569267272949218, .014715778641402721406, -.058503162115812301636, -.080628484487533569336, .013471723534166812896, .088432200253009796142, .0057333302684128284454, .053983848541975021362, .11002668738365173339, -.056655380874872207641, -.026702994480729103088, .093781538307666778564, -.054673958569765090942, -.13729652762413024902, -.018286351114511489868, .078907564282417297363, -.12712123990058898926, .029213106259703636169, -.040787089616060256958, .11390887945890426636, .067593149840831756591, -.071756012737751007081, -.079660005867481231689, -.040104631334543228149, .046582367271184921264, .090251527726650238037, .079574882984161376953, -.0071930391713976860046, .061825547367334365844, .14481072127819061279, -.080121636390686035156, .021467437967658042907, -.10324212908744812011, -.083766117691993713378, -.038257036358118057251, -.10646814107894897461, .18560406565666198731, .12912304699420928956, -.043887164443731307983, .049666363745927810668, .087872438132762908936, .090710788965225219726, .021344034001231193542, .035175830125808715821, .0058099604211747646331, -.039304725825786590576, -.030624080449342727661, .0016676177037879824638, -.15358072519302368164, -.0016108879353851079941, -.0061439997516572475433, -.10629564523696899414, -.048382952809333801269, .096326418220996856689, .046554282307624816894, -.00076445203740149736404, -.13558854162693023681, .076396353542804718017, -.12050648033618927001, .13775318861007690429, -.029819425195455551147, .022362027317285537719, -.050137817859649658203, .012289100326597690582, .088226906955242156982, -.00019194479682482779026, .064084306359291076661, .063600137829780578613, -.13246105611324310302, .017447149381041526794, .069785542786121368408, .061804030090570449829, -.041772555559873580932, -.060252070426940917968, .15938407182693481446, -.068411819636821746826, .097487136721611022949, .085471801459789276123, -.059593006968498229981, .043674904853105545043, -.20477126538753509521, .049175199121236801147, -.055517453700304031372, .10386756807565689086, -.11795367300510406494, .14392358064651489257, -.071758039295673370361, -.21455161273479461669, -.11355800926685333251, .012306759133934974671, .062935218214988708496, -.076733253896236419677, .080633088946342468261, -.11525595188140869141, .017615828663110733032, -.015210998244583606719, .041250526905059814453, -.013686502352356910706, -.35449263453483581542, .0026710352394729852676, .091578930616378784179, .20552384853363037109, -.021034220233559608459, .063262827694416046142, -.083335407078266143798, -.11970119923353195191, .012437067925930023193, .030780954286456108093, .022992219775915145874, .022368263453245162963, .066066578030586242676, -.084142744541168212891, .061110712587833404541, .20714576542377471923, .10172302275896072387, -.096266075968742370606, -.017330504953861236572, .078952215611934661866, -.068517580628395080566, -.0057096467353403568267, .10035095363855361938, .050398044288158416748, -.33835566043853759766, .044072706252336502076, -.023613989353179931641, .22618226706981658936, -.0090044457465410232543, -.14424835145473480224, -.037537563592195510864, -.012119695544242858886, .14918796718120574951, .078486897051334381103, .075012698769569396972, .0068048774264752864837, -.018572721630334854126, -.11508087068796157836, -.14797110855579376221, .10877022147178649902, -.038610812276601791381, -.085719779133796691894, .14319615066051483154, -.14610420167446136474, -.024826826527714729309, .16105084121227264404, .082391969859600067138, -.0084324143826961517333, -.046084817498922348022, .10822478681802749633, .056944396346807479858, -.12710206210613250732, .14380732178688049316, -.038927111774682998657, .060929212719202041626, .045901112258434295654, .052659131586551666259, .050329897552728652954, .055411703884601593017, .073711425065994262696, .072962552309036254882, -.0039118137210607528686, -.051354900002479553222, -.035267900675535202026, -.011948166415095329284, .075796641409397125244, -.0067796171642839908599, -.043093580752611160278, .000099209253676235675811, -.063098512589931488037, .061760183423757553101, -.013927742838859558106, .12463590502738952636, -.18979126214981079101, .024025280028581619262, .12755277752876281738, .018298316746950149536, .061236374080181121826, -.1117115020751953125, -.063949197530746459961, -.033463537693023681641, .057998713105916976928, -.15935491025447845458, .014615387655794620513, .055235195904970169067, .0034350070636719465256, -.10277443379163742066, .085209198296070098876, -.023268176242709159851, .085313640534877777099, .043945372104644775391, .035581987351179122924, .063343085348606109619, .015699539333581924438, -.019698953256011009216, .094785973429679870606, .15861411392688751221, -.21505700051784515381, -.061242200434207916259, .10139628499746322631, .054617173969745635986, .063631482422351837158, -.085129603743553161621, .12683668732643127441, -.070061668753623962402, .078650869429111480712, .19011779129505157471, -.066377580165863037109, .11143352836370468139, .20646293461322784423, .078123077750205993652, .10182949155569076538, .032863758504390716552, -.019082067534327507019, .12752631306648254394, .054198641330003738403, .12517069280147552491, -.057447381317615509033, -.040619760751724243164, -.057259462773799896241, .00086315307999029755592, -.053335383534431457519, .036726482212543487548, -.11905650794506072998, -.018268784508109092712, -.044493086636066436767, -.053184315562248229981, .019770724698901176452, .13655404746532440186, .086901746690273284912, .019807714968919754028, -.093643754720687866211, -.0065183024853467941284, .13562186062335968017, -.10695416480302810668, -.060802739113569259643, -.080440022051334381103, -.095137558877468109131, .15421338379383087158, .035559456795454025268, .044344145804643630981, -.13653728365898132324, .093051724135875701904, .14142389595508575439, -.061814371496438980102, -.060157086700201034546, -.078138031065464019776, .15519112348556518554, .045273046940565109252, -.061310347169637680053, -.075741097331047058106, -.13070251047611236572, -.14540983736515045166, -.11546280980110168457, -.093241296708583831787, .086008973419666290283, .0095527032390236854553, .031498748809099197387, -.052859056740999221801, -.11586606502532958984, -.062171354889869689941, .042789772152900695801, .059568256139755249023, -.15000912547111511231, -.058016795665025711059, -.033979550004005432128, .079871550202369689941, .15450774133205413818, .047726266086101531982, -.12199843674898147583, -.045893948525190353393, -.024430662393569946289, .066178776323795318603, -.059939287602901458741, -.091592229902744293212, -.18003596365451812744, .050483223050832748413, .14872919023036956787, -.12494672089815139771, -.037369720637798309326, .027677690610289573669, .070209480822086334228, -.088712714612483978271, .13036499917507171631, -.017092604190111160278, .025508837774395942687, .044629953801631927491, .10040614753961563111, -.077257245779037475586, .047057274729013442993, -.071298509836196899414, -.046474255621433258056, .16901618242263793946, .026769595220685005187, .086722590029239654541, -.066889874637126922607, .028233105316758155822, .12295199930667877197, -.16333736479282379151, .082589969038963317871, -.053739096969366073608, .092514663934707641601, .070871278643608093261, .11886804550886154174, -.053602315485477447509, .13626617193222045898, -.019748844206333160401, -.12089840322732925416, -.24070817232131958007, .15027897059917449951, -.038379933685064315796, -.074061460793018341064, -.15717884898185729981, -.036982282996177673339, -.014180591329932212829, .079939767718315124511, .17578063905239105224, -.10344955325126647949, .099240012466907501221, .13328579068183898926, .0050941850058734416961, .096172176301479339599, .17442426085472106933, -.057383131235837936401, -.041510514914989471436, -.030773533508181571961, -.060212571173906326293, -.042715799063444137573, -.0089184120297431945801, .057373613119125366211, .15370967984199523926, -.097737155854701995849, .076747842133045196533, .079159632325172424316, -.059677585959434509277, -.10163188725709915161, -.074042998254299163818, .051579102873802185058, .11258963495492935181, -.076367668807506561279, .077453128993511199951, -.042999591678380966186, -.039437878876924514771, -.025928571820259094238, .084689401090145111083, .073151111602783203125, .016415994614362716674, -.10574951022863388061, -.0088972672820091247558, .018638335168361663818, -.030070062726736068726, -.51254773139953613281, -.017395559698343276977, .060017820447683334351, -.048462286591529846191, -.070172399282455444336, -.039715271443128585816, .067449428141117095947, .069513797760009765625, -.044890139251947402954, -.11439596116542816162, .026105560362339019776, .094591103494167327881, -.060216549783945083618, .16226692497730255126, .022435313090682029724, -.11339483410120010376, .086511358618736267089, .14578899741172790527, -.074854061007499694824, .051716566085815429687, -.076930947601795196533, .045365013182163238526, -.034596573561429977416, -.072367243468761444091, .019366234540939331054, -.11108854413032531738, .060251396149396896362, -.12221951037645339966, -.058244105428457260131, .024633215740323066711, .085757680237293243408, .021583816036581993103, .0014467175351455807686, -.068062037229537963867, -.018024746328592300416, -.0078941294923424720764, -.045391503721475601196, -.073709778487682342529, -.072252906858921051026, .037061683833599090576, .13159501552581787109, -.055548656731843948364, .0026712866965681314468, -.066392727196216583251, .071197360754013061523, .073007613420486450196, -.058063946664333343506, .11829040199518203736, .090955629944801330566, .072909124195575714111, .068996280431747436523, .11048281192779541016, .085700117051601409912, -.073338046669960021972, -.049943100661039352416, .094952709972858428956, -.0038058375939726829528, .028341850265860557556, .10806329548358917236, -.021207498386502265931, .13731768727302551269, .13526929914951324462, .055056698620319366456, -.016316924244165420532, -.032298382371664047241, .022826373577117919921, .081599816679954528808, .015855750069022178649, .051428787410259246826, .012790031731128692626, .035054933279752731323, .0030938768759369850158, -.067092664539813995361, -.094230674207210540771, -.075678884983062744141, .27981230616569519042, -.080231681466102600097, -.10916028171777725219, -.031302411109209060668, .14119586348533630371, -.10244761407375335693, .11929044872522354126, .048006229102611541748, .073726765811443328857, .00041971082100644707679, -.062729284167289733886, -.19763784110546112061, .17552676796913146972, .12911927700042724609, .022614531219005584716, .051903288811445236206, -.074097551405429840087, .12041787803173065186, .12820836901664733886, .11165481805801391601, -.068390488624572753906, .00065528781851753592491, -.13809370994567871093, -.00022297918621916323901, -.084090411663055419921, .038129750639200210571, .015473307110369205474, -.13387905061244964599, -.11216911673545837402, .0090386644005775451661, .13227540254592895507, -.021621027961373329162, -.055052239447832107543, .11396015435457229614, -.068747542798519134521, -.063306257128715515136, -.26814174652099609375, .0037014798726886510848, -.017834832891821861267, .032523758709430694581, -.066172659397125244141, .037423323839902877807, .13851787149906158447, .10626149922609329223, -.0021064509637653827667, .0080012986436486244201, -.12753282487392425537, .026136932894587516784, .038572549819946289062, -.094934411346912384033, .026936430484056472778, .082451745867729187011, -.045962810516357421875, .10131878405809402466, .11180113255977630616, .11019238084554672241, .13204662501811981201, .018295768648386001586, .035303197801113128662, -.037810463458299636841, -.057661674916744232177, -.028056288138031959533, -.020935617387294769287, -.084718421101570129394, -.0031480379402637481689, -.062361415475606918334, -.26195785403251647949, -.034745510667562484741, .053431458771228790283, -.019200665876269340516, .14852149784564971923, .077810257673263549804, -.063412606716156005859, .088670611381530761718, .067361548542976379394, -.048671055585145950317, .089160487055778503417, -.039199512451887130737, .044708672910928726196, .022721663117408752441, .02366352081298828125, -.12433021515607833862, .064898967742919921875, -.010950516909360885621, -.011506100185215473176, .047743566334247589111, -.077690429985523223876, .053739793598651885986, .12521792948246002197, -.046849366277456283569, -.0089169908314943313598, -.086856290698051452636, .15187068283557891846, -.072726063430309295654, .019284209236502647399, -.037372548133134841918, .040343474596738815307, -.067755177617073059082, -.017583364620804786682, .019624615088105201721, .092085011303424835206, -.069585919380187988281, .021643171086907386779, -.080812290310859680176, .12826505303382873536, .13061548769474029541, .054610017687082290649, .061656333506107330322, -.19312471151351928711, -.038419298827648162841, .023252584040164947509, -.0069833658635616302491, -.17886374890804290771, .11288302391767501831, -.12335938215255737304, -.0072580138221383094787, .074872210621833801269, .10420162230730056762, -.044322747737169265747, .14676162600517272949, .0032918392680585384368, .14193791151046752929, -.014618578366935253143, -.00029767176602035760879, -.097034774720668792724, -.016640678048133850097, -.067487522959709167481, -.18384523689746856689, .018788164481520652771, -.059947635978460311889, -.041382052004337310791, -.064323648810386657714, -.038503907620906829833, -.0049685933627188205718, -.11444303393363952636, .015804948285222053527, -.083638072013854980468, .075691387057304382324, -.14261321723461151123, -.015377826057374477386, -.053159549832344055176, .044196400791406631469, -.076010003685951232911, -.048471234738826751708, -.067416213452816009521, -.051496837288141250611, -.018422512337565422058, -.10003858804702758789, -.015177638269960880279, .14476485550403594971, .12861829996109008789, .088524319231510162353, -.068521194159984588623, .12220853567123413086, .047840800136327743531, .013955649919807910919, -.034688089042901992797, -.072531796991825103759, -.011596234515309333801, .10722050815820693969, -.038290586322546005249, .056661039590835571289, -.076576113700866699218, .054596524685621261596, .067865602672100067138, .076995469629764556884, .23956349492073059082, -.11692409962415695191, -.0090113710612058639526, -.14294445514678955078, .036910176277160644531, .073884457349777221679, .043798502534627914428, -.047957174479961395263, -.052355796098709106446, -.039638582617044448852, .016472326591610908508, -.059484291821718215942, -.0023241108283400535583, -.056656558066606521606, .053520757704973220826, -.21285156905651092529, .051796186715364456176, .0067593785934150218963, -.073131494224071502686, .056344546377658843994, .032803203910589218139, -.045978084206581115722, .012059054337441921234, .074009522795677185058, .16344055533409118652, -.094489403069019317626, .089783981442451477051, .068552382290363311767, -.094720281660556793212, -.064518228173255920411, -.10681521147489547729, .13402818143367767333, .18673859536647796631, .15840637683868408203, -.10203147679567337036, -.10123059898614883422, .074032969772815704346, -.011930139735341072082, -.089628510177135467529, .038749642670154571533, .046832546591758728027, -.069618575274944305419, -.087661854922771453857, -.018450886011123657226, .045700557529926300048, .087684847414493560791, -.11261492222547531127, .061524435877799987792, -.10205468535423278808, -.039139948785305023193, .13543663918972015381, .059442967176437377929, -.14792498946189880371, .16928000748157501221, -.054070770740509033203, .13116887211799621582, .028384448960423469543, -.065740846097469329833, -.066516436636447906494, .030571749433875083923, .085575871169567108154, .069169595837593078613, .093178816139698028564, -.074192941188812255859, -.065236501395702362061, .052351437509059906006, .040969472378492355346, .0089074578136205673217, -.079830631613731384277, .031805172562599182128, .0047108801081776618957, -.16023847460746765136, .033315438777208328247, .20432311296463012696, .10753677040338516236, .054298993200063705444, .069054611027240753173, -.025415552780032157897, .031617492437362670898, -.12802569568157196044, -.13167838752269744873, .033167965710163116456, -.034123882651329040527, .014486155472695827484, .026442745700478553771, .054856471717357635498, .11334887892007827758, .13858899474143981933")
    # print(res)

    df = pd.read_csv('rows_for_db_final.csv', sep=';')
    # print(df.head())

    IrisDatabaseHandler().insert_embedded_sequence(df)

