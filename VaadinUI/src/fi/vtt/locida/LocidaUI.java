package fi.vtt.locida;

import com.vaadin.annotations.Theme;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.Sizeable;
import com.vaadin.server.VaadinRequest;
import com.vaadin.shared.ui.MarginInfo;
import com.vaadin.ui.Alignment;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Embedded;
import com.vaadin.ui.FormLayout;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Link;
import com.vaadin.ui.ListSelect;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.OptionGroup;
import com.vaadin.ui.Slider;
import com.vaadin.ui.TabSheet;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;
import com.vaadin.ui.themes.Reindeer;

import fi.vtt.RVaadin.RContainer;
import fi.vtt.RVaadin.RUpload;

/**
 * Main UI class
 */
@SuppressWarnings("serial")
@Theme("vtt")
public class LocidaUI extends UI {

	/* Application main window and the corresponding R instance */
	Window main;
	RContainer R;

	/* The analysis work flow is divided into separate tabs */
	private TabSheet maintabs = new TabSheet();
	private VerticalLayout submissionTab = new VerticalLayout();
	private VerticalLayout analysisTab = new VerticalLayout();
	private VerticalLayout downloadTab = new VerticalLayout();

	/* The lists to choose the treated (T) and control (C) variables from */
	private final ListSelect chooseT = new ListSelect("Treated (red)");
	private final ListSelect chooseC = new ListSelect("Control (blue)");

	/* Elements that need to be visible outside of their tab */
	private Button updateDetectionCurve = null;

	/* Which data set is selected for current analysis */
	private ComboBox dataSel;
	private final String[] existingDatasets = { "SR1" };

	@Override
	protected void init(VaadinRequest request) {
		final VerticalLayout main = new VerticalLayout();
		main.setMargin(true);
		setContent(main);

		R = new RContainer();

		R.eval("require(gdata)");
		R.eval("require(Locida)");

		/* Show PDF button for images */
		R.setGraphButtonsVisible(true);
		
		//
		// Build the User Interface
		//
		maintabs.setSizeFull();

		/* Data selection tab */
		buildSubmissionTab(submissionTab);
		maintabs.addTab(submissionTab, "Dataset");

		/* Analysis tab */
		buildAnalysisTab(analysisTab);
		maintabs.addTab(analysisTab, "Analysis");
		maintabs.getTab(analysisTab).setEnabled(false);

		/* Download tab */
		downloadTab = new VerticalLayout();
		Label download = new Label("Result Files");
		downloadTab.addComponent(download);
		maintabs.addTab(downloadTab, "Results");
		maintabs.getTab(downloadTab).setEnabled(false);

		main.addComponent(maintabs);

		/*
		 * Finally, close the application if the browser window is closed. Also
		 * clean up the R temporary files
		 */

		/* TODO: Make this work in Vaadin 7! */

		// main.addDetachListener(new Window.CloseListener() {
		//
		// public void windowClose(CloseEvent e) {
		// /* Close the R session */
		// try {
		// R.closeAndDeleteFiles();
		// } catch (Exception ex) {
		// ex.printStackTrace();
		// }
		//
		// /* Now close the Vaadin application window */
		// getMainWindow().getApplication().close();
		// }
		// });
	}

	/**
	 * Upload the Excel file and select the variables for analysis
	 * 
	 * @param submissionTab
	 */
	private void buildSubmissionTab(VerticalLayout submissionTab) {

		submissionTab.setMargin(true);
		submissionTab.setSizeFull();

		final TabSheet chooseOrSubmit = new TabSheet();
		chooseOrSubmit.setSizeFull();
		chooseOrSubmit.setStyleName(Reindeer.TABSHEET_MINIMAL);

		/* Read data from the database */
		HorizontalLayout dataBase = new HorizontalLayout();
		dataBase.setSpacing(true);
		dataBase.setMargin(true);

		dataSel = new ComboBox("Choose an Existing Dataset");

		for (String dataSetName : existingDatasets) {
			dataSel.addItem(dataSetName);
		}

		/* Put the first as default */
		dataSel.setValue(existingDatasets[0]);

		dataSel.setNullSelectionAllowed(false);
		dataBase.addComponent(dataSel);

		Button selectB = new Button("Select");
		dataBase.addComponent(selectB);
		dataBase.setComponentAlignment(selectB, Alignment.BOTTOM_RIGHT);

		chooseOrSubmit.addTab(dataBase, "Choose");

		/**
		 * Process the data from Database
		 */
		Button.ClickListener selectAction = new Button.ClickListener() {

			@Override
			public void buttonClick(ClickEvent event) {

				try {
					
					String dataSetName = dataSel.getValue().toString();
					
					R.tryEval("data('" + dataSetName + "')");
					R.tryEval("attach(" + dataSetName + ")");
					R.tryEval("nvd <- newLocida(X)");

					R.tryEval("factorVec <- samples[[factorName]]");
					R.tryEval("factors <- unique(factorVec)");

					/* Update the selection for treated and control groups */
					R.buildSelectOptions(chooseT, "factors");
					R.buildSelectOptions(chooseC, "factors");

					/* Enable and choose the analysis tab */
					maintabs.getTab(analysisTab).setEnabled(true);
					maintabs.setSelectedTab(analysisTab);

					/* Update the method preview image */
					updateDetectionCurve.click();

				} catch (Exception e) {
					e.printStackTrace();
					Notification.show("Cannot import the data set.");
				}
			}
		};
		selectB.addClickListener(selectAction);

		/* Upload data */
		VerticalLayout uploadPage = new VerticalLayout();
		uploadPage.setSpacing(true);
		uploadPage.setMargin(true);

		final RUpload uploadE = R.getUploadElement(null);
		uploadPage.addComponent(uploadE);

		final Button processB = new Button("Process the Selected File");

		uploadPage.addComponent(processB);
		chooseOrSubmit.addTab(uploadPage, "Submit File");

		uploadE.getSelection();

		/* Select the variables */
		final HorizontalLayout variablePage = new HorizontalLayout();
		variablePage.setSpacing(true);
		variablePage.setMargin(true);

		final ComboBox factorSel = new ComboBox("Divide into Factors Along");
		factorSel.setNullSelectionAllowed(false);

		final TextArea factorsOut = new TextArea("Factor Preview");
		factorsOut.setWidth("70ex");
		factorsOut.setWordwrap(false);
		final Button OKB = new Button("OK");

		variablePage.addComponent(factorSel);
		variablePage.addComponent(factorsOut);
		variablePage.addComponent(OKB);
		variablePage.setComponentAlignment(OKB, Alignment.BOTTOM_RIGHT);

		chooseOrSubmit.addTab(variablePage, "Set Factors");
		chooseOrSubmit.getTab(variablePage).setEnabled(false);

		submissionTab.addComponent(chooseOrSubmit);

		/**
		 * The actual upload processing
		 */
		Button.ClickListener processExcel = new Button.ClickListener() {

			public void buttonClick(ClickEvent event) {
				String selectedFile = uploadE.getSelection();

				/* The user didn't select any file */
				if (selectedFile == null) {
					Notification.show("Please choose one file from "
							+ "the \"Already submitted files\"");
					return;
				}

				/*
				 * Here, we let R to parse the submitted Excel file (using
				 * gdata). The data are stored in the statistical data.frame 'X'
				 * on the R side for further analysis.
				 */

				try {
					/* The R call is constructed by concatenating Strings */
					R.tryEval("W <- read.xls('" + selectedFile
							+ "', check.names=FALSE, stringsAsFactors=FALSE)");

					R.tryEval("suppressWarnings("
							+ "numericColNames <- as.numeric(colnames(W)))");

					/*
					 * This is a convention: Suppose that non-numeric column
					 * names correspond to factors, sample names and other
					 * descriptive data, and the rest is the spectrum data
					 */

					R.tryEval("factorColumns <- is.na(numericColNames);"
							+ "samples <- W[,factorColumns];"
							+ "X <- W[,!factorColumns]");

					/* Update the selection component to show the items */
					String[] colnames = R.tryEval("colnames(samples)")
							.asStrings();
					factorSel.removeAllItems();
					for (String colname : colnames) {
						factorSel.addItem(colname);
					}

					/*
					 * All went fine! We can already build the corresponding
					 * data structures in R, as it does depends only on X, not
					 * on the designFrame
					 */
					R.tryEval("nvd <- newLocida(X)");

					/* Let us enable and select the next tab */
					chooseOrSubmit.getTab(variablePage).setEnabled(true);
					chooseOrSubmit.setSelectedTab(variablePage);

					/* Update the method preview image */
					updateDetectionCurve.click();

				} catch (Exception e) {
					Notification.show("Cannot parse the Excel file.");
				}
			}
		};
		processB.addClickListener(processExcel);

		/**
		 * Listen the factor selection and update the Factor Preview
		 */
		ValueChangeListener factorSelected = new ValueChangeListener() {

			@Override
			public void valueChange(ValueChangeEvent event) {

				String factor = factorSel.getValue().toString();
				try {
					String[] factors = R.tryEval(
							"unique(samples[['" + factor + "']])").asStrings();

					StringBuffer sb = new StringBuffer();

					int i = 1;
					for (String s : factors) {
						sb.append(i++ + ": " + s + "\n");
					}
					factorsOut.setValue(sb.toString());

				} catch (Exception e) {
					/* How could we get here? */
					Notification.show("Cannot process the selected factor.");
				}

			}
		};
		factorSel.addValueChangeListener(factorSelected);
		factorSel.setImmediate(true);

		/**
		 * We have now chosen the factor name to divide the samples with. Let us
		 * tell that to R and the switch to the analysis tab.
		 */
		Button.ClickListener factorOK = new Button.ClickListener() {

			@Override
			public void buttonClick(ClickEvent event) {

				R.eval("factorName <- '" + factorSel.getValue().toString()
						+ "'");
				R.eval("factorVec <- samples[[factorName]]");
				R.eval("factors <- unique(factorVec)");

				/* Update the selection for treated and control groups */
				R.buildSelectOptions(chooseT, "factors");
				R.buildSelectOptions(chooseC, "factors");

				/*
				 * Build the corresponding R object /* Switch automatically to
				 * the analysis tab
				 */
				maintabs.getTab(analysisTab).setEnabled(true);
				maintabs.setSelectedTab(analysisTab);

			}
		};
		OKB.addClickListener(factorOK);
	}

	/**
	 * The actual data analysis
	 * 
	 * @param analysisTab
	 */
	private void buildAnalysisTab(VerticalLayout analysisTab) {

		analysisTab.setMargin(true);

		/*
		 * Construct an additional (sub)TabSheet to choose the factors, plot
		 * range, algorithm parameters, etc. The different analysis will be
		 * shown in a separate (sub)TabSheet
		 */

		TabSheet inputValueSheet = new TabSheet();
		inputValueSheet.setSizeFull();
		inputValueSheet.setStyleName(Reindeer.TABSHEET_MINIMAL);

		/* The factor selection elements */
		chooseT.setMultiSelect(true);
		R.buildSelectListener(chooseT, "Tfactors");
		chooseC.setMultiSelect(true);
		R.buildSelectListener(chooseC, "Cfactors");

		VerticalLayout groupsAndRange = new VerticalLayout();
		HorizontalLayout groups = new HorizontalLayout();
		groups.setSpacing(true);
		groups.setMargin(new MarginInfo(true, false, false, false));

		groups.addComponent(chooseT);
		groups.addComponent(chooseC);

		groupsAndRange.addComponent(groups);

		/* Controls for the graphics size and computation range */

		FormLayout graphicsForm = new FormLayout();
		graphicsForm.setCaption("Graphics Range and Size");
		graphicsForm.setSizeUndefined();

		R.eval("xlim <- c(10, 0.5)");
		R.eval("ylim <- c(0, 0.1)");
		R.eval("resolution <- c(800,500)");
		R.eval("ps <- 12");

		FormLayout range = new FormLayout();
		range.setSizeUndefined();

		HorizontalLayout xlim = R.getParameterLayout("Horizontal range (ppm):",
				"xlim");
		final HorizontalLayout ylim = R.getParameterLayout(
				"Vertical range (Intensity):", "ylim");
		HorizontalLayout res = R.getParameterLayout(
				"Image size (width x height):", "resolution");
		Slider ps = R.getSlider("Font Size:", "ps", 6, 20);
		ps.setWidth("10ex");

		range.addComponent(xlim);
		range.addComponent(ylim);
		range.addComponent(res);
		range.addComponent(ps);
		groupsAndRange.addComponent(range);

		inputValueSheet.addTab(groupsAndRange, "Grouping and Range");

		/* Detection parameter layout */

		final HorizontalLayout detectionPreview = new HorizontalLayout();

		R.eval("P <- c(0.05, 0.05)");
		R.eval("D <- c(0.001, 0.003)");

		R.eval("methods <- c('cut-off', 'linear', 'saturating')");
		R.eval("measures <- c('absolute', 'relative')");

		FormLayout detection = new FormLayout();
		detection.setSizeUndefined();
		detection.setMargin(new MarginInfo(false, true, false, false));

		OptionGroup measure = R.getOptionGroup("measures", "measure");
		measure.setCaption("Scale");
		measure.setValue("absolute");

		OptionGroup method = R.getOptionGroup("methods", "method");
		method.setCaption("Computation Method");
		method.setValue("saturating");

		detection.addComponent(measure);
		detection.addComponent(method);

		boolean[] enabledDefault = { true, false };
		final HorizontalLayout P = R.getParameterLayout("P-Value Threshold",
				"P", enabledDefault);
		final HorizontalLayout D = R.getParameterLayout("Detection Threshold",
				"D");
		detection.addComponent(P);
		detection.addComponent(D);

		/* Construct the preview window for the detection curve */
		final HorizontalLayout preview = new HorizontalLayout();
		final int previewSize = 300;
		preview.setWidth(previewSize, Sizeable.Unit.PIXELS);
		preview.setHeight(previewSize, Sizeable.Unit.PIXELS);

		/* A button to update the method preview graph */
		updateDetectionCurve = new Button("Update Preview",
				new Button.ClickListener() {

					@Override
					public void buttonClick(ClickEvent event) {
						/* Using no R closures would avoid this call */
						updateParameters();

						Embedded detectionCurve = R.getEmbeddedGraph(
								"nvd$getDetectionCurve()", previewSize,
								previewSize, "DetectionCurve");
						preview.removeAllComponents();
						preview.addComponent(detectionCurve);

					}
				});
		detection.addComponent(updateDetectionCurve);

		detectionPreview.addComponent(detection);
		detectionPreview.setComponentAlignment(detection,
				Alignment.MIDDLE_CENTER);

		detectionPreview.addComponent(preview);

		inputValueSheet.addTab(detectionPreview, "Detection Parameters");

		/**
		 * Add the buttons for parameter estimation
		 */

		/* Find the optimal Y range for the given x-axis */
		Button estimateYlim = new Button("Optimal Range",
				new Button.ClickListener() {

					@Override
					public void buttonClick(ClickEvent event) {
						if (updateDesign()) {
							R.eval("ylim <- nvd$getMaxYlim()");
							R.setParameterLayoutValues("ylim", ylim);

						}
					}
				});
		ylim.addComponent(estimateYlim);

		/* A checkbox to choose between a single and multiple p-values. */
		final CheckBox singlePThreshold = new CheckBox("Single Threshold", true);
		P.addComponent(singlePThreshold);

		ValueChangeListener pCheckBoxListener = new ValueChangeListener() {

			@Override
			public void valueChange(ValueChangeEvent event) {

				TextField[] p = R.getParameterLayoutTextFields(P);

				if ((Boolean) singlePThreshold.getValue()) {
					/* Set p1 same as p0 */
					p[1].setValue(R.getDouble("nvd$getP()[1]").toString());
					p[1].setEnabled(false);
				} else {
					/*
					 * Read the current p0, compute the Bonferroni-adjusted
					 * p-value and set that as p1
					 */
					R.eval("p0 <- nvd$getP()[1]");
					R.setParameterLayoutValues(
							R.getDoubles("c(p0, signif(p0/ncol(X),3))"), P);
					p[1].setEnabled(true);
				}
			}
		};
		singlePThreshold.addValueChangeListener(pCheckBoxListener);
		singlePThreshold.setImmediate(true);

		/* Trigger the D-value estimation */
		final Button estimateD = new Button("Estimate Thresholds",
				new Button.ClickListener() {

					@Override
					public void buttonClick(ClickEvent event) {
						if (updateDesign()) {
							R.eval("D <- nvd$getEstimateD()");
							R.setParameterLayoutValues("D", D);

							/* Update also the preview */
							updateDetectionCurve.click();
						}
					}
				});
		D.addComponent(estimateD);

		/*
		 * For measure and method, we automatically call the update by
		 * programmatically click():ing the update Button.
		 */
		ValueChangeListener methodChanged = new ValueChangeListener() {
			@Override
			public void valueChange(ValueChangeEvent event) {
				/* Check whether D can be re-estimated */
				if (!R.isSelectionEmpty(chooseT)
						&& !R.isSelectionEmpty(chooseC)) {
					estimateD.click();
				}
				/* In any case, the curve can be updated */
				updateDetectionCurve.click();
			}
		};
		method.addValueChangeListener(methodChanged);

		/* We also listen to measure and update the preview graph, if possible */
		ValueChangeListener measureChanged = new ValueChangeListener() {
			@Override
			public void valueChange(ValueChangeEvent event) {

				/* We update the detection curve only if D can be re-estimated */
				if (!R.isSelectionEmpty(chooseT)
						&& !R.isSelectionEmpty(chooseC)) {
					estimateD.click();
					updateDetectionCurve.click();

				}

			}
		};
		measure.addValueChangeListener(measureChanged);

		/**
		 * Layout and buttons for different analysis
		 */

		TabSheet runMethods = new TabSheet();
		runMethods.setStyleName(Reindeer.TABSHEET_MINIMAL);
		runMethods.setSizeFull();

		HorizontalLayout groupCompare = new HorizontalLayout();
		groupCompare.setMargin(true);
		groupCompare.setSpacing(true);

		Button computeGC = new Button("Compute");
		groupCompare.addComponent(computeGC);

		VerticalLayout export = new VerticalLayout();
		final Button exportGC = new Button("Export List");
		export.addComponent(exportGC);
		final CheckBox commaSep = new CheckBox("Decimal ,", false);
		export.addComponent(commaSep);
		groupCompare.addComponent(export);

		runMethods.addTab(groupCompare, "Group Comparison");

		VerticalLayout crossCompare = new VerticalLayout();
		crossCompare.setMargin(true);
		Button computeCC = new Button("Compute");
		crossCompare.addComponent(computeCC);

		runMethods.addTab(crossCompare, "Cross Comparison");

		analysisTab.addComponent(inputValueSheet);
		analysisTab.addComponent(runMethods);

		/**
		 * Triggers for the different analysis
		 */

		/* Group Comparison */
		Button.ClickListener drawGroupComparison = new Button.ClickListener() {

			@Override
			public void buttonClick(ClickEvent event) {

				if (updateDesign()) {
					Window plot = R.getGraph("par('ps'=ps); nvd$getPlot()",
							R.getInt("resolution[1]"),
							R.getInt("resolution[2]"));

					getUI().addWindow(plot);
				}
			}
		};
		computeGC.addClickListener(drawGroupComparison);

		Button.ClickListener exportListener = new Button.ClickListener() {
			@Override
			public void buttonClick(ClickEvent event) {

				if (updateDesign()) {

					/* Generate the result frame and construct the file name. */
					R.eval("rf <- nvd$getResultFrame();");
					R.eval("fileName <- nvd$getFileName()");

					/*
					 * The same, old problem with csv files and some Western
					 * European locales
					 */
					if ((Boolean) commaSep.getValue()) {
						R.eval("write.csv2(rf, fileName, row.names=FALSE)");
					} else {
						R.eval("write.csv(rf, fileName, row.names=FALSE)");
					}
					Link resultFile = R.getDownloadLink(
							R.getString("fileName"), R.getString("fileName"));

					maintabs.getTab(downloadTab).setEnabled(true);
					maintabs.setSelectedTab(downloadTab);
					downloadTab.addComponent(resultFile);
				}
			}
		};
		exportGC.addClickListener(exportListener);

		/* Cross-comparison */
		Button.ClickListener drawCrossComparison = new Button.ClickListener() {

			@Override
			public void buttonClick(ClickEvent event) {

				if (updateDesign()) {

					R.eval("res <- nvd$getDifferenceMatrices( factorVec, factors )");

					Window plot = R.getGraph("par('ps'=ps); "
							+ "nvd$getHeatmapDisplay( res, factors )",
							R.getInt("resolution[1]"),
							R.getInt("resolution[2]"));

					getUI().addWindow(plot);
				}
			}
		};
		computeCC.addClickListener(drawCrossComparison);
	}

	/**
	 * Because of the current analysis (nvd is a closure in R), we need to
	 * update it by hand. Make no closures, and you do not have to do this!
	 * 
	 * @return
	 */
	protected boolean updateParameters() {

		try {
			/* Update the R object */
			R.tryEval("nvd$setXlim( xlim );" + "nvd$setYlim( ylim );"
					+ "nvd$setP( P );" + "nvd$setD( D );"
					+ "nvd$setMeasure( measure );" + "nvd$setMethod( method )");

		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}

	/**
	 * This helper function updates the data structures beyond the simple global
	 * variable update made by the xVA User interface elements, including the
	 * construction of designFrame
	 * 
	 * @return boolean object
	 */
	protected boolean updateDesign() {

		/* Additional checks for the values and ranges */
		if (R.isSelectionEmpty(chooseT) || R.isSelectionEmpty(chooseC)) {
			Notification.show("Define the groups.", "Choose items from "
					+ "both of the lists.", Type.HUMANIZED_MESSAGE);
			return false;
		}
		/* Update the parameters */
		if (!updateParameters()) {
			return false;
		}

		try {
			/* Construct and set the designFrame */
			R.tryEval("designFrame <- data.frame("
					+ "'Treated' = samples[[factorName]] %in% Tfactors,"
					+ "'Control' = samples[[factorName]] %in% Cfactors)");
			R.tryEval("names(designFrame) <- c(paste(Tfactors, collapse='\n'),"
					+ "paste(Cfactors, collapse='\n'))");

			R.tryEval("nvd$setDesignFrame( designFrame )");

		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}
}
