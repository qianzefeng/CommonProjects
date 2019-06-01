#include "rtk_platform.h"
#include "rtk_trace.h"

/* execute command -------------------------------------------------------------
* execute command line by operating system shell
* args   : char   *cmd      I   command line
* return : execution status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int execcmd(const char *cmd)
{
#ifdef WIN32
	PROCESS_INFORMATION info;
	STARTUPINFO si={0};
	DWORD stat;
	char cmds[1024];

	trace(3,"execcmd: cmd=%s\n",cmd);

	si.cb=sizeof(si);
	sprintf(cmds,"cmd /c %s",cmd);
	if (!CreateProcess(NULL,(LPTSTR)cmds,NULL,NULL,FALSE,CREATE_NO_WINDOW,NULL,
		NULL,&si,&info)) return -1;
	WaitForSingleObject(info.hProcess,INFINITE);
	if (!GetExitCodeProcess(info.hProcess,&stat)) stat=-1;
	CloseHandle(info.hProcess);
	CloseHandle(info.hThread);
	return (int)stat;
#else
	trace(3,"execcmd: cmd=%s\n",cmd);

	return system(cmd);
#endif
}


/* expand file path ------------------------------------------------------------
* expand file path with wild-card (*) in file
* args   : char   *path     I   file path to expand (captal insensitive)
*          char   *paths    O   expanded file paths
*          int    nmax      I   max number of expanded file paths
* return : number of expanded file paths
* notes  : the order of expanded files is alphabetical order
*-----------------------------------------------------------------------------*/
extern int expath(const char *path, char *paths[], int nmax)
{
	int i,j,n=0;
	char tmp[1024];
#ifdef WIN32
	WIN32_FIND_DATA file;
	HANDLE h;
	char dir[1024]="",*p;

	trace(3,"expath  : path=%s nmax=%d\n",path,nmax);

	if ((p=strrchr(path,'\\'))) 
	{
		strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
	}

	if ((h=FindFirstFile((LPCTSTR)path,&file))==INVALID_HANDLE_VALUE) 
	{
		strcpy(paths[0],path);
		return 1;
	}
	sprintf(paths[n++],"%s%s",dir,file.cFileName);
	while (FindNextFile(h,&file)&&n<nmax) 
	{
		if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
		sprintf(paths[n++],"%s%s",dir,file.cFileName);
	}
	FindClose(h);
#else
	struct dirent *d;
	DIR *dp;
	const char *file=path;
	char dir[1024]="",s1[1024],s2[1024],*p,*q,*r;

	trace(3,"expath  : path=%s nmax=%d\n",path,nmax);

	if ((p=strrchr(path,'/'))||(p=strrchr(path,'\\'))) {
		file=p+1; strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
	}
	if (!(dp=opendir(*dir?dir:"."))) return 0;
	while ((d=readdir(dp))) {
		if (*(d->d_name)=='.') continue;
		sprintf(s1,"^%s$",d->d_name);
		sprintf(s2,"^%s$",file);
		for (p=s1;*p;p++) *p=(char)tolower((int)*p);
		for (p=s2;*p;p++) *p=(char)tolower((int)*p);

		for (p=s1,q=strtok_r(s2,"*",&r);q;q=strtok_r(NULL,"*",&r)) {
			if ((p=strstr(p,q))) p+=strlen(q); else break;
		}
		if (p&&n<nmax) sprintf(paths[n++],"%s%s",dir,d->d_name);
	}
	closedir(dp);
#endif
	/* sort paths in alphabetical order */
	for (i=0;i<n-1;i++) 
	{
		for (j=i+1;j<n;j++) 
		{
			if (strcmp(paths[i],paths[j])>0) 
			{
				strcpy(tmp,paths[i]);
				strcpy(paths[i],paths[j]);
				strcpy(paths[j],tmp);
			}
		}
	}
	for (i=0;i<n;i++) trace(3,"expath  : file=%s\n",paths[i]);

	return n;
}

/* create directory ------------------------------------------------------------
* create directory if not exist
* args   : char   *path     I   file path to be saved
* return : none
* notes  : not recursive. only one level
*-----------------------------------------------------------------------------*/
extern void createdir(const char *path)
{
	char buff[1024],*p;

	tracet(3,"createdir: path=%s\n",path);

	strcpy(buff,path);
	if (!(p=strrchr(buff,FILEPATHSEP))) return;
	*p='\0';

#ifdef WIN32
	CreateDirectoryA(buff,NULL);
#else
	mkdir(buff,0777);
#endif
}