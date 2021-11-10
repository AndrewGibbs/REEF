function [ExtraNodes, ExtraWeights, NodesToSkip] = QuadSqrtExtraPtNodes(order)


    if (order <= 1.5)
	M = [
	    1.172258571393266D-01  5.000000000000000D-01
	];
	NodesToSkip = 1;

    elseif (order <= 2)
	M = [
	    9.252112715421378D-02  4.198079625266162D-01
	    1.000000000000000D-00  1.080192037473384D+00
	];
	NodesToSkip = 2;

    elseif (order <= 2.5)
	M = [
	    6.023873796408450D-02  2.858439990420468D-01
	    8.780704050676215D-01  1.214156000957953D+00
	];
	NodesToSkip = 2;

    elseif (order <= 3)
	M = [
	    7.262978413470474D-03  3.907638767531813D-02
	    2.246325512521893D-01  4.873484056646474D-01
	    1.000000000000000D+00  9.735752066600344D-01
	];
	NodesToSkip = 2;

    elseif (order <= 3.5)
	M = [
	    1.282368909458828D-02  6.363996663105925D-02
	    2.694286346792474D-01  5.077434578043636D-01
	    1.018414523786358D+00  9.286165755645772D-01
	];
	NodesToSkip = 2;

    elseif (order <= 4)
	M = [
	    1.189242434021285D-02  5.927215035616424D-02
	    2.578220434738662D-01  4.955981740306228D-01
	    1.007750064585281D+00  9.427131290628058D-01
	    2.000000000000000D+00  1.002416546550407D+00
	];
	NodesToSkip = 3;

    elseif (order <= 6)
	M = [
	    3.317925942699451D-03  1.681780929883469D-02
	    8.283019705296352D-02  1.755244404544475D-01
	    4.136094925726231D-01  5.039350503858001D-01
	    1.088744373688402D+00  8.266241339680867D-01
	    2.006482101852379D+00  9.773065848981277D-01
	    3.000000000000000D+00  9.997919809947032D-01
	];
	NodesToSkip = 4;

    elseif (order <= 8)
	M = [
	    1.214130606523435D-03  6.199844884297793D-03
	    3.223952700027058D-02  7.106286791720044D-02
	    1.790935383649920D-01  2.408930104410471D-01
	    5.437663805244631D-01  4.975929263668960D-01
	    1.176116628396759D+00  7.592446540441226D-01
	    2.031848210716014D+00  9.322446399614420D-01
	    3.001961225690812D+00  9.928171438160095D-01
	    4.000000000000000D+00  9.999449125689846D-01
	];
	NodesToSkip = 5;

    elseif (order <= 10)
	M = [
	    1.745862989163252D-04  1.016950985948944D-03
	    8.613670540457314D-03  2.294670686517670D-02
	    6.733385088703690D-02  1.076657968022888D-01
	    2.514488774733840D-01  2.734577662465576D-01
	    6.341845573737690D-01  4.978815591924992D-01
	    1.248404055083152D+00  7.256208919565360D-01
	    2.065688031953401D+00  8.952638690320078D-01
	    3.009199358662542D+00  9.778157465381624D-01
	    4.000416269690208D+00  9.983390781399277D-01
	    5.000000000000000D+00  9.999916342408948D-01
	];
	NodesToSkip = 6;

    elseif (order <= 12)
	M = [
	    5.710218427206990D-04  2.921018926912141D-03
	    1.540424351115548D-02  3.431130611256885D-02
	    8.834248407196555D-02  1.224669495638615D-01
	    2.824462054509770D-01  2.761108242022520D-01
	    6.574869892305580D-01  4.797809643010337D-01
	    1.246541060977993D+00  6.966555677271379D-01
	    2.039218495130811D+00  8.790077941972658D-01
	    2.979333487049800D+00  9.868622449294327D-01
	    3.985772595393049D+00  1.015142389688201D+00
	    4.997240804311428D+00  1.006209712632210D+00
	    5.999868793951190D+00  1.000528829922287D+00
	    7.000000000000000D+00  1.000002397796838D+00
	];
	NodesToSkip = 8;

    elseif (order <= 14)
	M = [
	    3.419821460249725D-04  1.750957243202047D-03
	    9.296593430187960D-03  2.080726584287380D-02
	    5.406214771755252D-02  7.586830616433430D-02
	    1.763945096508648D-01  1.766020526671851D-01
	    4.218486605653738D-01  3.206624362072232D-01
	    8.274022895884040D-01  4.934405290553812D-01
	    1.410287585637014D+00  6.707497030698472D-01
	    2.160997505238153D+00  8.244959025366557D-01
	    3.043504749358223D+00  9.314646742162802D-01
	    4.005692579069439D+00  9.845768443163154D-01
	    4.999732707905968D+00  9.992852769154770D-01
	    5.999875191971098D+00  1.000273112957723D+00
	    6.999994560568667D+00  1.000022857402321D+00
	    8.000000000000000D+00  1.000000081405180D+00
	];
	NodesToSkip = 9;

    elseif (order <= 16)
	M = [
	    2.158438988280793D-04  1.105804873501181D-03
	    5.898432743709196D-03  1.324499944707956D-02
	    3.462795956896131D-02  4.899842307592144D-02
	    1.145586495070213D-01  1.165326192868815D-01
	    2.790344218856415D-01  2.178586693194957D-01
	    5.600113798653321D-01  3.481766016945031D-01
	    9.814091242883119D-01  4.964027915911545D-01
	    1.553594853974655D+00  6.469026189623831D-01
	    2.270179114036658D+00  7.823688971783889D-01
	    3.108234601715371D+00  8.877772445893361D-01
	    4.032930893996553D+00  9.551665077035583D-01
	    5.006803270228157D+00  9.876285579741800D-01
	    6.000815466735179D+00  9.979929183863017D-01
	    7.000045035079542D+00  9.998470620634641D-01
	    8.000000738923901D+00  9.999962891645340D-01
	    9.000000000000000D+00  9.999999946893169D-01
	];
	NodesToSkip = 10;
    else
	error(sprintf('QuadSqrtExtraPtNodes: order %d is too high', order));
    end

    ExtraNodes = M(:,1);
    ExtraWeights = M(:,2);

%end
