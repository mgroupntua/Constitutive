using MGroup.LinearAlgebra.Vectors;
using MGroup.Multiscale.SupportiveClasses;
using Xunit;

namespace MGroup.Multiscale.Tests
{
	public static class MicrostructureBvpNRNLAnalyzerTest
	{
		[Fact]
		public static void CheckMicrostructureBvpNRNLAnalyzer()
		{
			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MSmicroBvpAnalyzerTest\\U_sunol_micro_3.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MSmicroBvpAnalyzerTest\\U_sunol_micro_6.txt";
			double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
			double[] displacements2ndncrement = PrintUtilities.ReadVector(results_file2);

			(IVector uInitialFreeDOFs_state1, IVector uInitialFreeDOFs_state2) = NRNLAnalyzerDevelopTest.SolveDisplLoadsExample();


			Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1.CopyToArray()));
			Assert.True(NRNLAnalyzerDevelopTest.AreDisplacementsSame(displacements2ndncrement, uInitialFreeDOFs_state2.CopyToArray()));
		}
	}
}
