static double w[N] = {
  3.779211177528728e-18,
  5.035852722282574e-18,
  5.919542444491104e-18,
  6.625780997425304e-18,
  7.225554026797095e-18,
  7.753411501834412e-18,
  8.22901389351605e-18,
  8.66471612981438e-18,
  9.068837885372816e-18,
  9.447275165901862e-18,
  9.804374136660416e-18,
  1.0143440402100701e-17,
  1.046705311944029e-17,
  1.0777267798949829e-17,
  1.1075752268719609e-17,
  1.1363880762787455e-17,
  1.164280081091786e-17,
  1.1913481912009988e-17,
  1.2176751677985312e-17,
  1.2433323156881527e-17,
  1.2683815817413517e-17,
  1.2928771894948216e-17,
  1.3168669287123936e-17,
  1.3403931845028942e-17,
  1.3634937672194328e-17,
  1.386202588123086e-17,
  1.408550214314811e-17,
  1.4305643282033703e-17,
  1.452270110785588e-17,
  1.4736905636009493e-17,
  1.494846780931835e-17,
  1.5157581813408613e-17,
  1.5364427057490915e-17,
  1.5569169878084408e-17,
  1.57719650119731e-17,
  1.5972956875899105e-17,
  1.6172280683578905e-17,
  1.6370063425141637e-17,
  1.6566424729707162e-17,
  1.6761477628301348e-17,
  1.6955329231460542e-17,
  1.7148081333564528e-17,
  1.7339830954047733e-17,
  1.7530670824087173e-17,
  1.772068982608615e-17,
  1.7909973392213292e-17,
  1.8098603867375507e-17,
  1.8286660841268705e-17,
  1.8474221453535323e-17,
  1.8661360675541846e-17,
  1.8848151571856147e-17,
  1.9034665544139553e-17,
  1.922097255986149e-17,
  1.94071413679859e-17,
  1.959323970356166e-17,
  1.9779334482967482e-17,
  1.996549199141056e-17,
  2.015177806415392e-17,
  2.0338258262846282e-17,
  2.052499804824809e-17,
  2.071206295058628e-17,
  2.0899518738726243e-17,
  2.108743158932197e-17,
  2.127586825709298e-17,
  2.146489624737944e-17,
  2.1654583992144607e-17,
  2.1845001030626488e-17,
  2.203621819588908e-17,
  2.2228307808588562e-17,
  2.2421343879352476e-17,
  2.2615402321272053e-17,
  2.2810561174131224e-17,
  2.3006900842143298e-17,
  2.3204504347140732e-17,
  2.3403457599368684e-17,
  2.360384968827396e-17,
  2.3805773195962715e-17,
  2.4009324536330014e-17,
  2.421460432324996e-17,
  2.4421717771666483e-17,
  2.4630775135954096e-17,
  2.4841892190539375e-17,
  2.5055190758505437e-17,
  2.527079929476491e-17,
  2.548885353140859e-17,
  2.5709497194049845e-17,
  2.5932882799429966e-17,
  2.6159172546277296e-17,
  2.6388539313487046e-17,
  2.6621167782188563e-17,
  2.6857255701293365e-17,
  2.7097015319798532e-17,
  2.7340675013619677e-17,
  2.7588481140256304e-17,
  2.7840700161422706e-17,
  2.809762108226512e-17,
  2.835955826639796e-17,
  2.8626854699347176e-17,
  2.889988578991277e-17,
  2.917906382056743e-17,
  2.946484318580545e-17,
  2.9757726593421176e-17,
  3.0058272450907446e-17,
  3.0367103721563165e-17,
  3.068491861822208e-17,
  3.101250361502133e-17,
  3.1350749411368486e-17,
  3.170067069509885e-17,
  3.206343085064238e-17,
  3.2440373184106596e-17,
  3.283306085507585e-17,
  3.324332861780275e-17,
  3.3673350851217557e-17,
  3.412573248168003e-17,
  3.460363276607498e-17,
  3.511093738513032e-17,
  3.565250353172647e-17,
  3.623451883582939e-17,
  3.686504450680305e-17,
  3.75548699516312e-17,
  3.831892265142546e-17,
  3.9178734569779584e-17,
  4.0167091972264215e-17,
  4.1337714538177976e-17,
  4.2788395857263743e-17,
  4.472928950257001e-17,
  4.777594786313511e-17,
  5.1529423044407315e-17,
};

static uint64_t k[N] = {
  0,
  54076415620536904,
  61300587756530160,
  64377012539529264,
  66076294707096096,
  67151864522975000,
  67892968177924664,
  68434202989192568,
  68846593712197952,
  69171123658870632,
  69433082538507880,
  69648914306148520,
  69829769879067288,
  69983476194917440,
  70115687770302160,
  70230590878559368,
  70331350677449488,
  70420404420292040,
  70499659178058240,
  70570628466878520,
  70634528715794040,
  70692348697106704,
  70744900361805480,
  70792856639540568,
  70836779941489144,
  70877143927839008,
  70914350325727000,
  70948742062233352,
  70980613620876648,
  71010219282873288,
  71037779740400504,
  71063487444945400,
  71087510964150584,
  71109998555058080,
  71131081113284408,
  71150874621589976,
  71169482194156344,
  71186995792270104,
  71203497671331240,
  71219061606927376,
  71233753938254240,
  71247634459757496,
  71260757186039240,
  71273171010449816,
  71284920274100840,
  71296045259079904,
  71306582617266824,
  71316565744221304,
  71326025106040584,
  71334988525801840,
  71343481435146536,
  71351527095693856,
  71359146794246632,
  71366360015153064,
  71373184592683680,
  71379636845862160,
  71385731697832616,
  71391482781545128,
  71396902533285504,
  71402002275357128,
  71406792289034624,
  71411281878747928,
  71415479428314056,
  71419392449911592,
  71423027626383784,
  71426390847360208,
  71429487239599328,
  71432321191875536,
  71434896374660288,
  71437215754778184,
  71439281605151936,
  71441095509685592,
  71442658363269824,
  71443970366826712,
  71445031017241888,
  71445839091957312,
  71446392627917936,
  71446688894475496,
  71446724359753304,
  71446494649861640,
  71445994500222504,
  71445217698110208,
  71444157015336032,
  71442804129794104,
  71441149534335072,
  71439182431133664,
  71436890609354568,
  71434260303482624,
  71431276029149600,
  71427920392635256,
  71424173869414824,
  71420014546126464,
  71415417819088288,
  71410356040935128,
  71404798104978112,
  71398708954393280,
  71392049000153288,
  71384773427507192,
  71376831365482520,
  71368164886913096,
  71358707797302752,
  71348384158599168,
  71337106477510176,
  71324773465667520,
  71311267248289504,
  71296449855397552,
  71280158769685808,
  71262201219561744,
  71242346781854352,
  71220317675959912,
  71195775857061400,
  71168305596529440,
  71137389581305912,
  71102375511978928,
  71062428446027928,
  71016461188283864,
  70963029851113296,
  70900172224299768,
  70825148412611304,
  70734006393383296,
  70620815666934008,
  70476226522053056,
  70284534227628608,
  70017030195232984,
  69614581079045472,
  68930870409638744,
  67462502132127024,
  66808818195620008,
};

static double f[N] = {
  1.0,
  0.9635996931270862,
  0.9362826816850596,
  0.9130436479717402,
  0.8922816507840261,
  0.8732430489100695,
  0.8555006078694506,
  0.8387836052959896,
  0.822907211381409,
  0.8077382946829605,
  0.7931770117713051,
  0.7791460859296877,
  0.7655841738977045,
  0.7524415591746114,
  0.7396772436726473,
  0.7272569183441848,
  0.7151515074104986,
  0.7033360990161581,
  0.6917891434366751,
  0.6804918409973341,
  0.6694276673488904,
  0.658582000050088,
  0.6479418211102225,
  0.6374954773350423,
  0.6272324852499273,
  0.6171433708188809,
  0.6072195366251203,
  0.5974531509445167,
  0.5878370544347066,
  0.5783646811197631,
  0.5690299910679509,
  0.5598274127040869,
  0.5507517931146045,
  0.5417983550254255,
  0.5329626593838361,
  0.5242405726729841,
  0.5156282382440018,
  0.507122051075569,
  0.4987186354709795,
  0.4904148252838441,
  0.4822076463294852,
  0.47409430069301695,
  0.4660721526894561,
  0.45813871626787206,
  0.4502916436820392,
  0.44252871527546844,
  0.4348478302499909,
  0.4272469983049961,
  0.4197243320495744,
  0.412278040102661,
  0.40490642080722294,
  0.3976078564938733,
  0.3903808082373146,
  0.3832238110559012,
  0.3761354695105626,
  0.3691144536644722,
  0.3621594953693176,
  0.3552693848479171,
  0.3484429675463266,
  0.3416791412315504,
  0.3349768533135892,
  0.3283350983728503,
  0.3217529158759849,
  0.3152293880650109,
  0.3087636380061811,
  0.30235482778648354,
  0.296002156846933,
  0.28970486044295984,
  0.283462208223233,
  0.2772735029191881,
  0.2711380791383846,
  0.2650553022555892,
  0.25902456739620483,
  0.25304529850732577,
  0.2471169475123214,
  0.24123899354543982,
  0.23541094226347908,
  0.22963232523211613,
  0.22390269938500842,
  0.2182216465543054,
  0.2125887730717303,
  0.20700370943992652,
  0.20146611007431367,
  0.19597565311627774,
  0.19053204031913715,
  0.1851349970089922,
  0.17978427212329545,
  0.1744796383307895,
  0.169220892237365,
  0.16400785468342038,
  0.1588403711394793,
  0.15371831220818166,
  0.14864157424234226,
  0.14361008009062776,
  0.1386237799845946,
  0.13368265258343937,
  0.1287867061959432,
  0.12393598020286782,
  0.11913054670765083,
  0.11437051244886601,
  0.10965602101484027,
  0.10498725540942132,
  0.10036444102865587,
  0.09578784912173144,
  0.09125780082683026,
  0.08677467189478018,
  0.08233889824223566,
  0.0779509825139734,
  0.0736115018841134,
  0.06932111739357791,
  0.06508058521306807,
  0.060890770348040406,
  0.05675266348104985,
  0.05266740190305101,
  0.048636295859867805,
  0.044660862200491425,
  0.040742868074444175,
  0.0368843887866562,
  0.03308788614622575,
  0.02935631744000685,
  0.02569329193593427,
  0.022103304615927098,
  0.018592102737011288,
  0.015167298010546568,
  0.011839478657884862,
  0.008624484412859885,
  0.005548995220771345,
  0.002669629083880923,
};