function tests = testFIONAPopulationEFOScreen
    tests = functiontests(localfunctions);
end

function testAutomaticMedianAndStrictThreshold(testCase)
    screen = buildFIONAPopulationEFOScreen([90 95 100 100 105 150 151 NaN 0],1.5,'A',[]);
    verifyEqual(testCase,screen.ReferenceEFO,100);
    verifyEqual(testCase,screen.ReferenceTraceCount,7);
    verifyEqual(testCase,screen.ThresholdEFO,150);
    verifyEqual(testCase,find(screen.Flagged),7);
    verifyTrue(testCase,isnan(screen.Ratio(8)) && isnan(screen.Ratio(9)));
end

function testCapsPreserveScoreAndNeverPromote(testCase)
    trace = struct('Grade','S','Score',96,'Flags','none','NeedsManualReview',false);
    a = buildFIONAPopulationEFOScreen(160,1.5,'A',100);
    b = buildFIONAPopulationEFOScreen(160,1.5,'B',100);
    cappedA = applyFIONAPopulationEFOGradeCap(trace,a,1);
    cappedB = applyFIONAPopulationEFOGradeCap(trace,b,1);
    verifyEqual(testCase,cappedA.Grade,'A');
    verifyEqual(testCase,cappedB.Grade,'B');
    verifyEqual(testCase,cappedB.Score,96);
    verifyEqual(testCase,cappedB.PopulationEFO.UncappedGrade,'S');
    verifyTrue(testCase,cappedB.NeedsManualReview && cappedB.PopulationEFO.GradeCapped);
    verifyEqual(testCase,applyFIONAPopulationEFOGradeCap(cappedB,b,1),cappedB);
    trace.Grade = 'C';
    lower = applyFIONAPopulationEFOGradeCap(trace,b,1);
    verifyEqual(testCase,lower.Grade,'C');
    verifyFalse(testCase,lower.PopulationEFO.GradeCapped);
    verifyTrue(testCase,lower.PopulationEFO.Flagged);
end

function testDisabledAndInsufficientData(testCase)
    off = buildFIONAPopulationEFOScreen([100 500],1.5,'off',100);
    small = buildFIONAPopulationEFOScreen([100 500],1.5,'A',[]);
    empty = buildFIONAPopulationEFOScreen(zeros(0,1),1.5,'A',[]);
    verifyFalse(testCase,any(off.Flagged));
    verifyFalse(testCase,any(small.Flagged));
    verifyEmpty(testCase,empty.Flagged);
    verifyEqual(testCase,small.ReferenceSource,'insufficient eligible traces');
end
