# Streamlit Cloud Deployment Guide

## Prerequisites

1. **Repository Structure**: Ensure all data files are included in the repository
2. **File Size Limits**: Be aware of Streamlit Cloud's file size restrictions
3. **Memory Constraints**: Monitor memory usage during large analyses

## Required Files

### Data Files (must be in repository):
```
data/
├── backgrounds/
│   ├── all_genes.txt
│   └── alias.json
├── libraries/
│   ├── *.gmt files
│   └── alias.json
├── gene_lists/
│   └── example_gene_list.txt
└── recent_release/
    └── gene_history

Homo_sapiens.gene_info
```

### Configuration Files:
- `requirements.txt` - Python dependencies
- `packages.txt` - System dependencies
- `.streamlit/config.toml` - Streamlit configuration

## Deployment Steps

1. **Push to GitHub**: Ensure all files are committed and pushed
2. **Connect to Streamlit Cloud**: Link your GitHub repository
3. **Configure Settings**:
   - Python version: 3.12
   - Main file path: `code/streamlit_app.py`
4. **Deploy**: Click deploy and monitor for errors

## Potential Issues & Solutions

### Memory Issues
- **Symptom**: App crashes during large analyses
- **Solution**: Implement data streaming or chunking for large datasets

### File System Issues
- **Symptom**: Cannot write results files
- **Solution**: Use Streamlit's session state for temporary storage

### Timeout Issues
- **Symptom**: Long-running analyses fail
- **Solution**: Implement progress indicators and break down large tasks

## Performance Optimization

1. **Lazy Loading**: Load data files only when needed
2. **Caching**: Use `@st.cache_data` for expensive computations
3. **Progress Indicators**: Show progress for long-running operations
4. **Error Handling**: Implement robust error handling for cloud environment

## Monitoring

- Check Streamlit Cloud logs for errors
- Monitor memory usage during peak times
- Test with various input sizes
- Verify all features work in cloud environment

