/*
 * Find synapes from Vglut dots Homer 
 * Author Philippe Mailly
 */

import fiji.util.gui.GenericDialogPlus;
import ij.*;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import mcib3d.geom.Objects3DPopulationColocalisation;
import mcib3d.geom.PairColocalisation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import Vglut_Stardist.StarDist2D;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.FilenameFilter;
import mcib3d.image3d.ImageLabeller;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;

public class Vglut_Homer_2D implements PlugIn {

    private boolean canceled = false;
    private String imageDir = "";
    private String outDirResults = "";
    private Calibration cal = new Calibration();
    
    private BufferedWriter outPutResults;
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public double minDots = 0.01;
    public double maxDots = 0.5;
    public String[] dotsDetectors = {"StarDist","DOG"};
    public String dotsDetector = "StarDist";
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    
    // Stardist
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThreshDots = 0.7;
    public final double stardistOverlayThreshDots = 0.35;
    public String stardistOutput = "Label Image"; 
    public String starDistModel = "";
    
    // DOG detections
    private String[] thMeths = AutoThresholder.getMethods();
    private String ch0ThMet = thMeths[10];
    public int ch0DogRadius1 = 2;
    public int ch0DogRadius2 = 4;
    private String ch1ThMet = thMeths[3];
    public int ch1DogRadius1 = 1;
    public int ch1DogRadius2 = 2;
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    
   
    
    
    /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }


    private boolean dialog() {
        canceled = true;
        String[] models = findStardistModels();
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 180, 0);
        gd.addImage(icon);
        gd.addDirectoryField("Choose images directory", imageDir);
        gd.addMessage("Dots detection method", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Dots segmentation method :",dotsDetectors, dotsDetectors[0]);
        gd.addMessage("StarDist model", Font.getFont("Monospace"), Color.blue);
        if (models.length >= 2) {
            gd.addChoice("StarDist model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist model :", starDistModel);
        }
        gd.addMessage("DOG dots detection", Font.getFont("Monospace"), Color.blue);
        gd.addMessage("Channel 0", Font.getFont("Monospace"), Color.black);
        gd.addChoice("Threshold method   :",thMeths, ch0ThMet);
        gd.addNumericField("DOG radius 1 :", ch0DogRadius1);
        gd.addNumericField("DOG radius 2 :", ch0DogRadius2);
        gd.addMessage("Channel 1", Font.getFont("Monospace"), Color.black);
        gd.addChoice("Threshold method   :",thMeths, ch1ThMet);
        gd.addNumericField("DOG radius 1 :", ch1DogRadius1);
        gd.addNumericField("DOG radius 2 :", ch1DogRadius2);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = false;
        imageDir = gd.getNextString()+ File.separator;
        dotsDetector = gd.getNextChoice();
        if (models.length >= 2) {
            starDistModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            starDistModel = gd.getNextString();
        }
        if (dotsDetector.equals("StarDist") && starDistModel.isEmpty()) {
            IJ.error("No model specify !!");
            canceled = true;
        }
        ch0ThMet = gd.getNextChoice();
        ch0DogRadius1 = (int)gd.getNextNumber();
        ch0DogRadius2 = (int)gd.getNextNumber();
        ch1ThMet = gd.getNextChoice();
        ch1DogRadius1 = (int)gd.getNextNumber();
        ch1DogRadius2 = (int)gd.getNextNumber();
        return(canceled);
    }
    
    
    /*
    Find starDist models in Fiji models folder
    */
    private String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        return(models);
    }
    
    
     /**  
     * median 3D box filter
     * Using CLIJ2
     * @param imgCL
     * @param sizeX
     * @param sizeY
     * @param sizeZ
     * @return imgOut
     */ 
    public ClearCLBuffer medianFilter(ClearCLBuffer imgCL, double sizeX, double sizeY, double sizeZ) {
        ClearCLBuffer imgIn = clij2.push(imgCL);
        ClearCLBuffer imgOut = clij2.create(imgIn);
        clij2.median2DBox(imgIn, imgOut, sizeX, sizeY);
        clij2.release(imgCL);
        return(imgOut);
    }
    
     /**
     * Find dots population with Stardist
     */
    public Objects3DPopulation stardistDotsPop(ImagePlus imgDots) throws IOException{
        ImagePlus img = new Duplicator().run(imgDots);
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLMed = medianFilter(imgCL, 1, 1, 1);
        clij2.release(imgCL);
        ImagePlus imgDotsMed = clij2.pull(imgCLMed);
        clij2.release(imgCLMed);
        // Go StarDist
        File starDistModelFile = new File(starDistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgDotsMed); 
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThreshDots, stardistOverlayThreshDots, stardistOutput);
        star.run();
        Img<? extends RealType<?>> img1 = star.label.getImgPlus().getImg();
        ImagePlus imgLab = ImageJFunctions.wrap((RandomAccessibleInterval)img1, "Labelled");
        imgLab.setCalibration(cal);
        Objects3DPopulation dotsPop = new Objects3DPopulation(imgLab);
        ArrayList<Object3D> objectsWithinVolume = dotsPop.getObjectsWithinVolume(minDots, maxDots, true);
        Objects3DPopulation dotsFilter = new Objects3DPopulation(objectsWithinVolume);
        imgDotsMed.close();
        imgLab.close();
        return(dotsFilter);
    } 
    
    /** 
     * Find dots with DOG CLIJ Method
     * @param img channel
     * @return dots population
     */
    public Objects3DPopulation findDotsDoGCLIJ(ImagePlus img, int dogRadius1, int dogRadius2, String thMet) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = DOG(imgCL, dogRadius1, dogRadius2);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(threshold(imgCLDOG, thMet));
        clij2.release(imgCLDOG);
        imgBin.setCalibration(img.getCalibration());
        Objects3DPopulation pmlPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minDots, maxDots, true));
        imgBin.close();
        return(pmlPop);
    } 
    
    
    public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
   
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        return(imgCLBin);
    }
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
        
     /**
     * Find synapses 
     * if vglut dots and homer dots have one commun pixel
     * else remove vglut and homer object
     * keep only objects that make synapses
     * @param vglutPop
     * @param homerPop
     * @return 
     */
    
    private ArrayList<Point3D> findSynapses (Objects3DPopulation vglutDots, Objects3DPopulation homerDots) {
        Objects3DPopulationColocalisation coloc = new Objects3DPopulationColocalisation(vglutDots,homerDots);
        ArrayList<PairColocalisation> pairColoc = coloc.getAllColocalisationPairs();
        ArrayList<Point3D> synapsesPop = new ArrayList();
        for (PairColocalisation pair : pairColoc) {
            if (pair.getVolumeColoc() > 0) {
                Object3D vglutObj = pair.getObject3D1();
                vglutObj.setValue(1);
                Object3D homerObj = pair.getObject3D2();
                homerObj.setValue(1);
                Point3D synapse = new Point3D((vglutObj.getCenterX() + homerObj.getCenterX())/2, (vglutObj.getCenterY() + homerObj.getCenterY())/2, 1);
                synapsesPop.add(synapse);
            }
        }
        return(synapsesPop);
    }
    
    

    
    /**
     * Draw dots
     */
    private void drawDots(ArrayList<Point3D> dots, ImagePlus img){
        for (int i = 0; i < dots.size(); i++) {
            int r = 6;
            int x = dots.get(i).getRoundX() - r/2;
            int y = dots.get(i).getRoundY() - r/2;
            int z = dots.get(i).getRoundZ() + 1;
            img.setSlice(z);
            ImageProcessor ip = img.getProcessor();
            ip.setColor(255);
            ip.drawOval(x, y, r, r);
            img.updateAndDraw();
        }
    }
    
    /**
     * Save objects Population in image
     * @param vglutPop
     * @param homerPop
     * @param synapsePop
     * @param imageName
     * @param img 
     */
    private void saveImgObjects(Objects3DPopulation vglutPop, Objects3DPopulation homerPop, ArrayList<Point3D> synapsePoints, String imageName, ImagePlus img) {
        //create image objects population
        ImageHandler imgVglutObj = ImageInt.wrap(img).createSameDimensions();
        ImageHandler imgHomerObj = imgVglutObj.duplicate();
        ImagePlus imgSynObj = imgHomerObj.duplicate().getImagePlus();
        
        //vglut population green
        if (vglutPop.getNbObjects() > 0)
            vglutPop.draw(imgVglutObj, 255);
        //homer population red
        if (homerPop.getNbObjects() > 0)
            homerPop.draw(imgHomerObj, 255);
        //Synapses population blue      
        if (synapsePoints.size() > 0)
            drawDots(synapsePoints, imgSynObj);
        
        //save image
        ImagePlus[] imgColors = {imgHomerObj.getImagePlus(), imgVglutObj.getImagePlus(), imgSynObj};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imageName+".tif");
    }
    
    /**
     * Find mean Area dots
    */
    private DescriptiveStatistics findArea(Objects3DPopulation pop) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < pop.getNbObjects(); i++) {
            Object3D obj = pop.getObject(i);
            if (obj.getValue() == 1)
                stats.addValue(obj.getAreaUnit());
        }
        return(stats);
    }
    
    
    /** Write results
     * 
     */
    private void writeResults(String image, Objects3DPopulation vglutPop, Objects3DPopulation homerPop, ArrayList<Point3D> synapses) throws IOException {
        DescriptiveStatistics vglutStats = findArea(vglutPop);
        DescriptiveStatistics homerStats = findArea(homerPop);
        outPutResults.write(image+"\t"+vglutPop.getNbObjects()+"\t"+vglutStats.getMean()+"\t"+vglutStats.getStandardDeviation()+"\t"+homerPop.getNbObjects()+
                "\t"+homerStats.getMean()+"\t"+homerStats.getStandardDeviation()+"\t"+synapses.size()+"\n");
        outPutResults.flush();
    }
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        FileWriter fileResults = null;
        try {
            if (!checkInstalledModules()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            if (!dialog()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            fileResults = new FileWriter(outDirResults + "results.xls", false);
            Arrays.sort(imageFile);
            int imageNum = 0;
            for (String f : imageFile) {
                // open image files
                if (f.endsWith(".msr")) {
                    imageNum++;
                    try {
                        // for all image files
                        // get file name without extension
                        String rootName = f.replace(".msr", "");
                        if (imageNum == 1) {
                            outPutResults = new BufferedWriter(fileResults);
                            outPutResults.write("Image Name\tVglut number\tVglut Mean Area\tVglut Std Area\tHomer number\tHomer Mean Area\tHomer Std Area"
                                    + "\tSynapse number\n");
                            outPutResults.flush();
                        }
                        // Find spatial calibration
                        // create OME-XML metadata store of the latest schema version
                        ServiceFactory factory;
                        factory = new ServiceFactory();
                        OMEXMLService service = factory.getInstance(OMEXMLService.class);
                        IMetadata meta = service.createOMEXMLMetadata();
                        ImageProcessorReader reader = new ImageProcessorReader();
                        reader.setMetadataStore(meta);
                        String imageName = imageDir+f;
                        reader.setId(imageName);
                        reader.setSeries(0);
                        System.out.println("Image "+rootName);
                        
                        /* read channels
                        *  series0 = Vglut
                        *  series1 = Homer
                        */
                        // Open channels
                        // Find homer
                        
                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        // Find homer
                        options.setSeriesOn(1, true);
                        ImagePlus imgHomer = BF.openImagePlus(options)[0];
                        cal = imgHomer.getCalibration();
                        Objects3DPopulation homerPop = new Objects3DPopulation();
                        if (dotsDetector.equals("StarDist")) {
                            IJ.showStatus("Finding homer dots with StarDist...");
                            homerPop = stardistDotsPop(imgHomer);
                        }
                        else {
                            IJ.showStatus("Finding homer dots with DOG ...");
                            homerPop = findDotsDoGCLIJ(imgHomer, ch1DogRadius1, ch1DogRadius2, ch1ThMet);
                        }
                        System.out.println(homerPop.getNbObjects()+" homer dots");
                        
                        // Find vglut
                        options.setSeriesOn(0, true);
                        ImagePlus imgVglut = BF.openImagePlus(options)[0];
                        Objects3DPopulation vglutPop = new Objects3DPopulation();
                        if (dotsDetector.equals("StarDist")) {
                            IJ.showStatus("Finding Vglut dots with StarDist...");
                            vglutPop = stardistDotsPop(imgVglut);
                        }
                        else {
                            IJ.showStatus("Finding Vglut dots with DOG...");
                            vglutPop = findDotsDoGCLIJ(imgVglut, ch1DogRadius1, ch0DogRadius2, ch0ThMet);
                        }
                        System.out.println(vglutPop.getNbObjects()+" vglut dots");
                        // Find synapses
                        IJ.showStatus("Finding synapses ...");
                        ArrayList<Point3D> synapses = findSynapses(vglutPop, homerPop);
                        System.out.println(synapses.size()+" synapses found");
                        
                        // Save image objects
                        saveImgObjects(vglutPop, homerPop, synapses, outDirResults+rootName, imgHomer);
                        
                        imgVglut.close();
                        imgHomer.close();
                        
                        writeResults(rootName, vglutPop, homerPop, synapses);
                        
                    } catch (DependencyException | ServiceException | FormatException | IOException ex) {
                        Logger.getLogger(Vglut_Homer_2D.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
            try {
                outPutResults.close();
            } catch (IOException ex) {
                Logger.getLogger(Vglut_Homer_2D.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showStatus("Process done");
        } catch (IOException ex) {
            Logger.getLogger(Vglut_Homer_2D.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }    
}
