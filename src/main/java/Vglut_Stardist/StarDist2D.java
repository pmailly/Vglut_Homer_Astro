package Vglut_Stardist;


import static Vglut_Stardist.Utils.raiToDataset;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

import javax.swing.JOptionPane;

import org.scijava.command.Command;
import org.scijava.command.CommandModule;

//import de.csbdresden.CommandFromMacro;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;


public class StarDist2D extends StarDist2DBase implements Command {
    private Dataset input;
    private boolean normalizeInput = true;
    private double percentileBottom = 0.2;
    private double percentileTop = 99.8;
    private Dataset prob;
    private Dataset dist;
    public Dataset label;
    private double probThresh = 0.55;
    private double nmsThresh = 0.4;
    private String outputType = "Image Label";
   
    private int nTiles = 1;
    private int excludeBoundary = 2;  // boundary_exclusion
    private String roiPosition = "Automatic";
    private String roiPositionActive = null;
    private boolean verbose = false;
    private boolean showCsbdeepProgress = false;
    private boolean showProbAndDist = false;
    private ImageJ ij;
    private Object obj_;
    private File starDistModelFile = null;

   
    
    private void checkForCSBDeep() {
        try {
            Class.forName("de.csbdresden.csbdeep.commands.GenericNetwork");
        } catch (ClassNotFoundException e) {
            JOptionPane.showMessageDialog(null,
                    "<html><p>"
                    + "StarDist relies on the CSBDeep plugin for neural network prediction.<br><br>"
                    + "Please install CSBDeep by enabling its update site.<br>"
                    + "Go to <code>Help > Update...</code>, then click on <code>Manage update sites</code>.<br>"
                    + "Please see <a href='https://github.com/csbdeep/csbdeep_website/wiki/CSBDeep-in-Fiji-%E2%80%93-Installation'>https://tinyurl.com/csbdeep-install</a> for more details."
                    + "</p><img src='"+getResource("images/csbdeep_updatesite.png")+"' width='498' height='324'>"
                    ,
                    "Required CSBDeep plugin missing",
                    JOptionPane.ERROR_MESSAGE);
            throw new RuntimeException("CSBDeep not installed");
        }
    }
    

    public StarDist2D(Object obj, File starDistModel) {
         ij = new ImageJ();
         ij.launch();
        dataset = ij.dataset();
        command = ij.command();
        obj_ = obj;
        starDistModelFile = starDistModel;
    }
    
    
    /**
     * Check image size
     * change nTiles if image > 2048x2048 
     */
    
    public void checkImgSize(ImagePlus img) {
        int width = img.getWidth(); 
        if (img.getWidth() >= 2048)
            nTiles = Math.round(width / 2048);
        else {
            nTiles = 1;
        }
    }

    

    @Override
    public void run() {
        checkForCSBDeep();
        if (!checkInputs()) return;

        if (roiPosition.equals("Automatic"))
            roiPositionActive = input.numDimensions() > 3 && !input.isRGBMerged() ? "Hyperstack" : "Stack";
        else
            roiPositionActive = roiPosition;

        try {
            final HashMap<String, Object> paramsCNN = new HashMap<>();
            paramsCNN.put("input", input);
            paramsCNN.put("normalizeInput", true);
            paramsCNN.put("percentileBottom", 0.2);
            paramsCNN.put("percentileTop", 99.8);
            paramsCNN.put("clip", false);
            paramsCNN.put("nTiles", nTiles);
            paramsCNN.put("blockMultiple", 64);
            paramsCNN.put("overlap", 64);
            paramsCNN.put("batchSize", 1);
            paramsCNN.put("showProgressDialog", false);
            paramsCNN.put("modelFile", starDistModelFile); 
            

            final HashMap<String, Object> paramsNMS = new HashMap<>();
            paramsNMS.put("probThresh", 0.55);
            paramsNMS.put("nmsThresh", 0.4);
            paramsNMS.put("excludeBoundary", 2);
            paramsNMS.put("roiPosition", roiPositionActive);
            paramsNMS.put("verbose", false);      
            
            final LinkedHashSet<AxisType> inputAxes = Utils.orderedAxesSet(input);
            final boolean isTimelapse = inputAxes.contains(Axes.TIME);

            // TODO: option to normalize image/timelapse channel by channel or all channels jointly
            
            if (true && isTimelapse) {
                // TODO: option to normalize timelapse frame by frame (currently) or jointly
                final ImgPlus<? extends RealType<?>> inputImgPlus = input.getImgPlus();
                final long numFrames = input.getFrames();
                final int inputTimeDim = IntStream.range(0, inputAxes.size()).filter(d -> input.axis(d).type() == Axes.TIME).findFirst().getAsInt();
                for (int t = 0; t < numFrames; t++) {
                    final Dataset inputFrameDS = Utils.raiToDataset(dataset, "Input Frame",
                            Views.hyperSlice(inputImgPlus, inputTimeDim, t),
                            inputAxes.stream().filter(axis -> axis != Axes.TIME));
                    paramsCNN.put("input", inputFrameDS);
                    final Future<CommandModule> futureCNN = command.run(de.csbdresden.csbdeep.commands.GenericNetwork.class, false, paramsCNN);
                    final Dataset prediction = (Dataset) futureCNN.get().getOutput("output");

                    final Pair<Dataset, Dataset> probAndDist = splitPrediction(prediction);
                    final Dataset probDS = probAndDist.getA();
                    final Dataset distDS = probAndDist.getB();
                    paramsNMS.put("prob", probDS);
                    paramsNMS.put("dist", distDS);
                    paramsNMS.put("outputType", outputType);
                    if (showProbAndDist) {
                        if (t==0) log.error(String.format("\"%s\" not implemented/supported for timelapse data.", "Show CNN Output"));
                    }

                    final Future<CommandModule> futureNMS = command.run(StarDist2DNMS.class, false, paramsNMS);
                    final Candidates polygons = (Candidates) futureNMS.get().getOutput("polygons");
                    export(outputType, polygons, 1+t, numFrames, roiPositionActive);

                    IJ.showProgress(1+t, (int)numFrames);
                }
                label = labelImageToDataset(outputType);                
                // if (roiManager != null) OverlayCommands.listRois(roiManager.getRoisAsArray());

            } else {
                // note: the code below supports timelapse data too. differences to above:
                //       - joint normalization of all frames
                //       - requires more memory to store intermediate results (prob and dist) of all frames
                //       - allows showing prob and dist easily
                final Future<CommandModule> futureCNN = command.run(de.csbdresden.csbdeep.commands.GenericNetwork.class, false, paramsCNN);
                final Dataset prediction = (Dataset) futureCNN.get().getOutput("output");

                final Pair<Dataset, Dataset> probAndDist = splitPrediction(prediction);
                final Dataset probDS = probAndDist.getA();
                final Dataset distDS = probAndDist.getB();
                paramsNMS.put("prob", probDS);
                paramsNMS.put("dist", distDS);
                paramsNMS.put("outputType", outputType);
                if (showProbAndDist) {
                    prob = probDS;
                    dist = distDS;
                }

                final Future<CommandModule> futureNMS = command.run(StarDist2DNMS.class, false, paramsNMS);
                label = (Dataset) futureNMS.get().getOutput("label");
            }
            // call at the end of the run() method
            //CommandFromMacro.record(this, this.command);
            
        } catch (InterruptedException | ExecutionException  e) {
            e.printStackTrace();
        }
    }

    // this function is very cumbersome... is there a better way to do this?
    private Pair<Dataset, Dataset> splitPrediction(final Dataset prediction) {
        final RandomAccessibleInterval<FloatType> predictionRAI = (RandomAccessibleInterval<FloatType>) prediction.getImgPlus();
        final LinkedHashSet<AxisType> predAxes = Utils.orderedAxesSet(prediction);

        final int predChannelDim = IntStream.range(0, predAxes.size()).filter(d -> prediction.axis(d).type() == Axes.CHANNEL).findFirst().getAsInt();
        final long[] predStart = predAxes.stream().mapToLong(axis -> {
            return axis == Axes.CHANNEL ? 1 : 0;
        }).toArray();
        final long[] predSize = predAxes.stream().mapToLong(axis -> {
            return axis == Axes.CHANNEL ? prediction.dimension(axis)-1 : prediction.dimension(axis);
        }).toArray();

        final RandomAccessibleInterval<FloatType> probRAI = Views.hyperSlice(predictionRAI, predChannelDim, 0);
        final RandomAccessibleInterval<FloatType> distRAI = Views.offsetInterval(predictionRAI, predStart, predSize);

        final Dataset probDS = Utils.raiToDataset(dataset, "Probability/Score Image", probRAI, predAxes.stream().filter(axis -> axis != Axes.CHANNEL));
        final Dataset distDS = Utils.raiToDataset(dataset, "Distance Image", distRAI, predAxes);

        return new ValuePair<>(probDS, distDS);
    }


    private boolean checkInputs() {
        final Set<AxisType> axes = Utils.orderedAxesSet(input);
        if (!( (input.numDimensions() == 2 && axes.containsAll(Arrays.asList(Axes.X, Axes.Y))) ||
               (input.numDimensions() == 3 && axes.containsAll(Arrays.asList(Axes.X, Axes.Y, Axes.TIME))) ||
               (input.numDimensions() == 3 && axes.containsAll(Arrays.asList(Axes.X, Axes.Y, Axes.CHANNEL))) ||
               (input.numDimensions() == 4 && axes.containsAll(Arrays.asList(Axes.X, Axes.Y, Axes.CHANNEL, Axes.TIME))) ))
            return showError("Input must be a 2D image or timelapse (with or without channels).");
        
        return true;
    }


    @Override
    protected void exportPolygons(Candidates polygons) {}


    @Override
    protected ImagePlus createLabelImage() {
        return IJ.createImage(outputType, "16-bit black", (int)input.getWidth(), (int)input.getHeight(), 1, 1, (int)input.getFrames());
    }


    public static void main(final String... args) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.launch(args);

        Dataset input = ij.scifio().datasetIO().open(StarDist2D.class.getClassLoader().getResource("yeast_crop.tif").getFile());
//        Dataset input = ij.scifio().datasetIO().open(StarDist2D.class.getClassLoader().getResource("yeast_timelapse.tif").getFile());
//        Dataset input = ij.scifio().datasetIO().open(StarDist2D.class.getClassLoader().getResource("patho_hyperstack.tif").getFile());
        ij.ui().show(input);
        
//        Recorder recorder = new Recorder();
//        recorder.show();

        final HashMap<String, Object> params = new HashMap<>();
        ij.command().run(StarDist2D.class, true, params);
    }
    
    public void loadInput(ImagePlus imp) {
        String modeldir = IJ.getDirectory("imagej")+"/models/";
        if ( imp.getNSlices()>1) imp.setDimensions(1, 1, imp.getNSlices());
        final AxisType[] axes = new AxisType[]{Axes.X, Axes.Y, Axes.TIME};
        final Img inputImg = (Img) ImageJFunctions.wrap(imp);
        input = raiToDataset(dataset, "input", inputImg, axes);      
    }
    
    public void setParams(double percentileBottomVar, double percentileTopVar, double probThreshVar, double overlapThreshVar, String outPutType){
        percentileBottom = percentileBottomVar;
        percentileTop = percentileTopVar;
        probThresh = probThreshVar;
        nmsThresh = overlapThreshVar;
        outputType = outPutType;
    }

}
